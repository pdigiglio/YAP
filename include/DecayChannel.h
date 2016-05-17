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

#ifndef yap_DecayChannel_h
#define yap_DecayChannel_h

#include "fwd/DecayChannel.h"

#include "fwd/CachedDataValue.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayingParticle.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/Model.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/ParticleCombination.h"
#include "fwd/SpinAmplitude.h"
#include "fwd/StatusManager.h"

#include "AmplitudeComponent.h"
#include "ReportsParticleCombinations.h"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \class DecayChannel
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald
class DecayChannel :
    public AmplitudeComponent,
    public ReportsParticleCombinations
{
public:

    /// \typedef map_type
    /// \brief maps SpinAmplitude-shared_ptr -> (map of spin projection -> free amplitude)
    using map_type = SpinAmplitudeMap<std::map<int, std::shared_ptr<ComplexParameter> > >;

    /// \name Constructors
    /// @{

    /// N-particle Constructor (at the moment only valid for 2 particles).
    /// DecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    DecayChannel(const ParticleVector& daughters);

    /// @}

    /// Calculate complex amplitude
    /// \param d DataPoint to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param two_m 2 * the spin projection to calculate for
    /// \param sm StatusManager to update
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                           int two_m, StatusManager& sm) const;

    /// check consistency of object
    virtual bool consistent() const override;

    /// \return vector of shared_ptr's to final-state particles of channel (recursively checked)
    std::vector<std::shared_ptr<FinalStateParticle> > finalStateParticles() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const ParticleVector& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude objects
    SpinAmplitudeVector spinAmplitudes();

    /// Get SpinAmplitude objects (const)
    const SpinAmplitudeVector spinAmplitudes() const
    { return const_cast<DecayChannel*>(this)->spinAmplitudes(); }

    /// Get map of (spin projection)->(free amplitude) for spin amplitude
    map_type::mapped_type& amplitudes(const std::shared_ptr<SpinAmplitude>& sa)
    { return Amplitudes_.at(sa); }

    /// Get map of (spin projection)->(free amplitude) for spin amplitude (const)
    const map_type::mapped_type& amplitudes(const std::shared_ptr<SpinAmplitude>& sa) const
    { return Amplitudes_.at(sa); }

    /// \return vector of ParticleCombinations
    const ParticleCombinationVector& particleCombinations() const
    { return ParticleCombinations_; }

    /// \return raw pointer to model through first Daughter
    const Model* model() const;

    /// \return raw pointer to owning DecayingParticle
    DecayingParticle* decayingParticle() const
    { return DecayingParticle_; }

    /// @}

    /// add a spin amplitude
    void addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa);

    /// \return free amplitude
    /// \param two_M twice the spin projection
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin
    std::shared_ptr<ComplexParameter> freeAmplitude(int two_M, unsigned l, unsigned two_s);

    /// \return Vector of free amplitudes
    ComplexParameterVector freeAmplitudes();

    /// Grant friend status to DecayingParticle to set itself as owner
    /// and to call fixSolitaryFreeAmplitudes()
    friend DecayingParticle;

protected:

    /// Add particle combination
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> c);

    /// call fixSolitaryFreeAmplitudes() on each daughter
    void fixSolitaryFreeAmplitudes();

    /// set raw pointer to owning DecayingParticle
    void setDecayingParticle(DecayingParticle* dp);

private:

    /// daughters of the decay
    ParticleVector Daughters_;

    /// Map of SpinAmplitude (by shared_ptr) to map(spin projection -> free amplitude)
    map_type Amplitudes_;

    /// raw pointer owning DecayingParticle
    DecayingParticle* DecayingParticle_;

    /// vector of shared_ptr<ParticleCombination>
    ParticleCombinationVector ParticleCombinations_;

};

/// convert to string
std::string to_string(const DecayChannel& dc);

/// << operator
inline std::ostream& operator<< (std::ostream& os, const DecayChannel& dc)
{ os << to_string(dc); return os; }

}

#endif
