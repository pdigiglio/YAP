/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/// \file

#ifndef yap_UnitSpinAmplitude_h
#define yap_UnitSpinAmplitude_h

#include "fwd/Spin.h"

#include "Constants.h"
#include "SpinAmplitude.h"

namespace yap {

/// \class UnitSpinAmplitude
/// \brief Implements a spin amplitude that always returns Complex_1
/// \author Daniel Greenwald
/// \ingroup SpinAmplitude
/// \todo Remove actual caching
class UnitSpinAmplitude : public SpinAmplitude
{
public:

    /// \return Complex_1 regardless of input
    const std::complex<double> calc(int two_M, const SpinProjectionVector& two_m,
                                    const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return Complex_1; }

    /// \return "unit-valued"
    virtual std::string formalism() const
    { return "unit-valued"; }

    /// grant SpinAmplitudeCache friend status to call constructor
    friend class SpinAmplitudeCache;

protected:

    /// Constructor
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    UnitSpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s)
        : SpinAmplitude(two_J, two_j, l, two_s)
    {
        addAmplitude(0, SpinProjectionVector(two_j.size(), 0));
    }

};

}

#endif