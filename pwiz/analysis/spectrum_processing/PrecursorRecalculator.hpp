//
// PrecursorRecalculator.hpp
//
//
// Original author: Darren Kessner <Darren.Kessner@cshs.org>
//
// Copyright 2008 Spielberg Family Center for Applied Proteomics
//   Cedars-Sinai Medical Center, Los Angeles, California  90048
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


#ifndef _PRECURSORRECALCULATOR_HPP_ 
#define _PRECURSORRECALCULATOR_HPP_ 


#include "data/msdata/MSData.hpp"


namespace pwiz {
namespace analysis {


class PWIZ_API_DECL PrecursorRecalculator
{
    public:

    struct PrecursorInfo
    {
        double mz;
        double intensity;
        double charge;
        double score;

        PrecursorInfo() : mz(0), intensity(0), charge(0), score(0) {}
    };

    virtual void recalculate(const pwiz::msdata::MZIntensityPair* begin,
                             const pwiz::msdata::MZIntensityPair* end,
                             const PrecursorInfo& initialEstimate,
                             std::vector<PrecursorInfo>& result) = 0;

    virtual ~PrecursorRecalculator() {}
}; 


} // namespace analysis 
} // namespace pwiz


#endif // _PRECURSORRECALCULATOR_HPP_ 

