//
// $Id$
//
//
// Original author: Matt Chambers <matt.chambers .@. vanderbilt.edu>
//
// Copyright 2014 Vanderbilt University - Nashville, TN 37232
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//


#define PWIZ_SOURCE

#pragma unmanaged
#include "ShimadzuReader.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/DateTime.hpp"


#pragma managed
#include "pwiz/utility/misc/cpp_cli_utilities.hpp"
using namespace pwiz::util;


using System::String;
using System::Math;
using System::Collections::Generic::IList;
//namespace ShimadzuAPI = Shimadzu::LabSolutions::DataReader;
namespace ShimadzuIO = Shimadzu::LabSolutions::IO;
namespace ShimadzuGeneric = ShimadzuIO::Generic;
using ShimadzuIO::Data::DataObject;
using ShimadzuIO::Method::MethodObject;
typedef ShimadzuGeneric::Tool ShimadzuUtil;
//using ShimadzuAPI::ReaderResult;

namespace pwiz {
namespace vendor_api {
namespace Shimadzu {


class ChromatogramImpl : public Chromatogram
{
    public:
    ChromatogramImpl(ShimadzuIO::Generic::MassChromatogramObject^ chromatogram, const SRMTransition& transition)
        : transition_(transition),
          chromatogram_(chromatogram)
    {}

    virtual const SRMTransition& getTransition() const { return transition_; }
    virtual int getTotalDataPoints() const { try { return (int) chromatogram_->TotalPoints; } CATCH_AND_FORWARD }
    virtual void getXArray(std::vector<double>& x) const
    {
        try
        {
            auto timeArray = chromatogram_->RetTimeList;
            x.resize(timeArray->Length);
            for (size_t i = 0; i < x.size(); ++i)
                x[i] = timeArray[i] / 1000.0;
        } CATCH_AND_FORWARD
    }

    virtual void getYArray(std::vector<double>& y) const
    {
        try
        {
            auto intensityArray = chromatogram_->ChromIntList;
            y.resize(intensityArray->Length);
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = intensityArray[i];
        } CATCH_AND_FORWARD
    }

    private:
    SRMTransition transition_;
    gcroot<ShimadzuIO::Generic::MassChromatogramObject^> chromatogram_;
};


class SpectrumImpl : public Spectrum
{
public:
    SpectrumImpl(ShimadzuIO::Generic::MassSpectrumObject^ spectrum)
        : spectrum_(spectrum)
    {}

    virtual double getScanTime() const { return spectrum_->RetentionTime; }
    virtual int getMSLevel() const { return spectrum_->PrecursorMzList->Count == 0 ? 1 : 2; }

    virtual bool getHasIsolationInfo() const { return false; }
    virtual void getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const { }

    virtual bool getHasPrecursorInfo() const { return spectrum_->PrecursorMzList->Count > 0; }
    virtual void getPrecursorInfo(double& selectedMz, double& intensity, int& charge) const
    {
        if (!getHasPrecursorInfo())
            return;

        selectedMz = spectrum_->PrecursorMzList[0] / ShimadzuUtil::MASSNUMBER_UNIT;
        intensity = 0;
        charge = spectrum_->PrecursorChargeState;
    }

    virtual int getTotalDataPoints(bool doCentroid) const { try { return doCentroid ? spectrum_->CentroidList->Count : spectrum_->ProfileList->Count; } CATCH_AND_FORWARD }
    virtual void getProfileArrays(std::vector<double>& x, std::vector<double>& y) const
    {
        try
        {
            auto profileArray = spectrum_->ProfileList;
            x.resize(profileArray->Count);
            y.resize(x.size());
            for (size_t i = 0; i < x.size(); ++i)
            {
                auto point = profileArray[i];
                x[i] = point->Mass / ShimadzuUtil::MASSNUMBER_UNIT;
                y[i] = point->Intensity;
            }
        } CATCH_AND_FORWARD
    }

    virtual void getCentroidArrays(std::vector<double>& x, std::vector<double>& y) const
    {
        try
        {
            auto centroidArray = spectrum_->CentroidList;
            x.resize(centroidArray->Count);
            y.resize(x.size());
            for (size_t i = 0; i < x.size(); ++i)
            {
                auto point = centroidArray[i];
                x[i] = point->Mass / ShimadzuUtil::MASSNUMBER_UNIT;
                y[i] = point->Intensity;
            }
        } CATCH_AND_FORWARD
    }

    private:
    gcroot<ShimadzuIO::Generic::MassSpectrumObject^> spectrum_;
};


class ShimadzuReaderImpl : public ShimadzuReader
{
    public:
    ShimadzuReaderImpl(const string& filepath)
    {
        try
        {
            dataObject_ = gcnew DataObject();
            String^ systemFilepath = ToSystemString(filepath);
            auto result = dataObject_->IO->LoadData(systemFilepath);
            if (ShimadzuUtil::Failed(result))
                throw runtime_error("[ShimadzuReader::ctor] LoadData error: " + ToStdString(System::Enum::GetName(result.GetType(), (System::Object^) result)));

            /*methodObject_ = gcnew MethodObject();
            result = methodObject_->IO->LoadMethod(systemFilepath);
            if (ShimadzuUtil::Failed(result))
                throw runtime_error("[ShimadzuReader::ctor] LoadMethod error: " + ToStdString(System::Enum::GetName(result.GetType(), (System::Object^) result)));*/

            auto chromatogramMng = dataObject_->MS->Chromatogram;

            scanCount_ = 0;
            segmentCount_ = chromatogramMng->SegmentCount;
            eventNumbersBySegment_.resize(segmentCount_);
            for (int i = 1; i <= segmentCount_; ++i)
            {
                auto& eventNumbers = eventNumbersBySegment_[i-1];
                eventNumbers.resize(chromatogramMng->EventCount(i));
                for (int j = 1; j <= eventNumbers.size(); ++j)
                {
                    eventNumbers[j - 1] = chromatogramMng->GetEventNo(i, j);

                    auto eventTIC = gcnew ShimadzuGeneric::MassChromatogramObject();
                    chromatogramMng->GetTICChromatogram(eventTIC, i, j);
                    scanCount_ += eventTIC->TotalPoints;
                }
            }
        }
        CATCH_AND_FORWARD
    }

    virtual ~ShimadzuReaderImpl() { dataObject_->IO->Close(); }

    virtual int getScanCount() const { return scanCount_; }

    //virtual std::string getVersion() const = 0;
    //virtual DeviceType getDeviceType() const = 0;
    //virtual std::string getDeviceName(DeviceType deviceType) const = 0;
    //virtual boost::local_time::local_date_time getAcquisitionTime() const = 0;

    virtual const set<SRMTransition>& getTransitions() const
    {
        if (!transitionSet_.empty())
            return transitionSet_;

        try
        {
            System::Collections::Generic::List<ShimadzuGeneric::Param::MS::MassEventInfo^>^ eventList;
            auto result = dataObject_->MS->Parameters->GetEventInfo(eventList);
            //auto result = methodObject_->MS->Parameters->GetEventInfo(eventList);
            //if (ShimadzuUtil::Failed(result))
            //    throw runtime_error("[ShimadzuReader::getTransitions] GetEventInfo error: " + ToStdString(System::Enum::GetName(result.GetType(), (System::Object^) result)));

            for each (auto evt in eventList)
            {
                if (evt->AnalysisMode != ShimadzuGeneric::AcqModes::MRM)
                    continue;
                SRMTransition t;
                t.channel = evt->Channel;
                t.event = evt->Event;
                t.segment = evt->Segment;
                t.collisionEnergy = evt->CE; // always non-negative, even if scan polarity is negative
                t.polarity = (short) evt->Polarity;
                t.Q1 = evt->StartMz / 10000.0;
                t.Q3 = evt->EndMz / 10000.0;
                transitionSet_.insert(transitionSet_.end(), t);
                transitions_[make_pair(t.segment, t.event)] = evt;
            }
            return transitionSet_;
        }
        CATCH_AND_FORWARD
    }

    virtual boost::local_time::local_date_time getAnalysisDate(bool adjustToHostTime) const
    {
        System::DateTime acquisitionTime = dataObject_->SampleInfo->AnalysisDate;
        bpt::ptime pt(boost::gregorian::date(acquisitionTime.Year, boost::gregorian::greg_month(acquisitionTime.Month), acquisitionTime.Day),
            bpt::time_duration(acquisitionTime.Hour, acquisitionTime.Minute, acquisitionTime.Second, bpt::millisec(acquisitionTime.Millisecond).fractional_seconds()));

        if (adjustToHostTime)
        {
            bpt::time_duration tzOffset = bpt::second_clock::universal_time() - bpt::second_clock::local_time();
            return blt::local_date_time(pt + tzOffset, blt::time_zone_ptr()); // treat time as if it came from host's time zone; actual time zone may not be provided by DateTime
        }
        else
            return blt::local_date_time(pt, blt::time_zone_ptr());
    }

    virtual ChromatogramPtr getChromatogram(const SRMTransition& transition) const
    {
        auto t = gcnew ShimadzuGeneric::MzTransition();
        t->Segment = transition.segment-1;
        t->Event = transition.event;
        t->Channel = transition.channel-1;
        t->StartMass = transition.Q1 * 10000;
        t->EndMass = transition.Q3 * 10000;
        try
        {
            ShimadzuGeneric::MassChromatogramObject^ chromatogram_;
            auto result = dataObject_->MS->Chromatogram->GetChromatogrambyEvent(chromatogram_, t);
            if (chromatogram_->TotalPoints == 0 && ShimadzuUtil::Failed(result))
                throw runtime_error("[ShimadzuReader::getChromatogram] GetChromatogrambyEvent error: " + ToStdString(System::Enum::GetName(result.GetType(), (System::Object^) result)));
            return ChromatogramPtr(new ChromatogramImpl(chromatogram_, transition));
        } CATCH_AND_FORWARD
    }

    virtual SpectrumPtr getSpectrum(int scanNumber) const
    {
        ShimadzuGeneric::MassSpectrumObject^ spectrum;
        dataObject_->MS->Spectrum->GetMSSpectrumByScan(spectrum, scanNumber);
        return SpectrumPtr(new SpectrumImpl(spectrum));
    }

    private:
    gcroot<DataObject^> dataObject_;
    //gcroot<MethodObject^> methodObject_;
    //gcroot<ShimadzuGeneric::MassChromatogramObject^> tic_;
    int segmentCount_;
    int scanCount_;
    vector<vector<int>> eventNumbersBySegment_;
    mutable map<pair<short, short>, gcroot<ShimadzuGeneric::Param::MS::MassEventInfo^> > transitions_;
    mutable set<SRMTransition> transitionSet_;
};


PWIZ_API_DECL
ShimadzuReaderPtr ShimadzuReader::create(const string& filepath)
{
    try { return ShimadzuReaderPtr(new ShimadzuReaderImpl(filepath)); } CATCH_AND_FORWARD
}


#pragma unmanaged
PWIZ_API_DECL
bool SRMTransition::operator< (const SRMTransition& rhs) const
{
    return (channel == rhs.channel ? ((segment == rhs.segment) ? (event < rhs.event) : (segment < rhs.segment)) : channel < rhs.channel);
}


} // Shimadzu
} // vendor_api
} // pwiz

