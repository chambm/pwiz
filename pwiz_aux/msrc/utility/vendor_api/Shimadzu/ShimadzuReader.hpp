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


#ifndef _SHIMADZUREADER_HPP_
#define _SHIMADZUREADER_HPP_


#include "pwiz/utility/misc/Export.hpp"
#include <string>
#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <boost/date_time.hpp>


namespace pwiz {
namespace vendor_api {
namespace Shimadzu {


struct PWIZ_API_DECL TimeRange { double start, end; };


struct PWIZ_API_DECL SRMTransition
{
    short channel;
    short event;
    short segment;
    double collisionEnergy;
    short polarity;
    double Q1;
    double Q3;
    //TimeRange acquiredTimeRange;

    bool operator< (const SRMTransition& rhs) const;
};


struct PWIZ_API_DECL Chromatogram
{
    virtual const SRMTransition& getTransition() const = 0;
    virtual int getTotalDataPoints() const = 0;
    virtual void getXArray(std::vector<double>& x) const = 0;
    virtual void getYArray(std::vector<double>& y) const = 0;

    virtual ~Chromatogram() {}
};

typedef boost::shared_ptr<Chromatogram> ChromatogramPtr;


struct PWIZ_API_DECL Spectrum
{
    virtual double getScanTime() const = 0;
    virtual int getMSLevel() const = 0;

    virtual bool getHasIsolationInfo() const = 0;
    virtual void getIsolationInfo(double& centerMz, double& lowerLimit, double& upperLimit) const = 0;

    virtual bool getHasPrecursorInfo() const = 0;
    virtual void getPrecursorInfo(double& selectedMz, double& intensity, int& charge) const = 0;

    virtual int getTotalDataPoints(bool doCentroid = false) const = 0;
    virtual void getProfileArrays(std::vector<double>& x, std::vector<double>& y) const = 0;
    virtual void getCentroidArrays(std::vector<double>& x, std::vector<double>& y) const = 0;

    virtual ~Spectrum() {}
};

typedef boost::shared_ptr<Spectrum> SpectrumPtr;


class PWIZ_API_DECL ShimadzuReader
{
public:
    typedef boost::shared_ptr<ShimadzuReader> Ptr;
    static Ptr create(const std::string& filepath);

    //virtual std::string getVersion() const = 0;
    //virtual DeviceType getDeviceType() const = 0;
    //virtual std::string getDeviceName(DeviceType deviceType) const = 0;
    virtual boost::local_time::local_date_time getAnalysisDate(bool adjustToHostTime) const = 0;

    virtual const std::set<SRMTransition>& getTransitions() const = 0;
    virtual ChromatogramPtr getChromatogram(const SRMTransition& transition) const = 0;

    //virtual ChromatogramPtr getTIC() const = 0;

    virtual int getScanCount() const = 0;
    virtual SpectrumPtr getSpectrum(int scanNumber) const = 0;

    virtual ~ShimadzuReader() {}
};

typedef ShimadzuReader::Ptr ShimadzuReaderPtr;


} // Shimadzu
} // vendor_api
} // pwiz


#endif // _SHIMADZUREADER_HPP_
