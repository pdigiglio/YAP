#include "Resonance.h"

#include "DecayChannel.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::unique_ptr<MassShape>& massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(nullptr)
{
    MassShape_.swap(massShape);
}

//-------------------------
std::complex<double> Resonance::amplitude(DataPartition& d, std::shared_ptr<const ParticleCombination> pc) const
{
    return DecayingParticle::amplitude(d, pc) * MassShape_->amplitude(d, pc);
}

//-------------------------
bool Resonance::consistent() const
{
    bool consistent = true;

    consistent &= DecayingParticle::consistent();
    consistent &= MassShape_->consistent();

    if (! consistent) {
        LOG(ERROR) << "Resonance is not consistent:  " << this->name() << "\n";
    }

    return consistent;
}

//-------------------------
void Resonance::setInitialStateParticle(InitialStateParticle* isp)
{
    DecayingParticle::setInitialStateParticle(isp);
    // hand ISP to mass shape
    MassShape_->setInitialStateParticle(initialStateParticle());
}

//-------------------------
void Resonance::addChannel(std::unique_ptr<DecayChannel>& c)
{
    for (auto& pc : c->particleCombinations())
        MassShape_->addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(pc));

    DecayingParticle::addChannel(c);
}

//-------------------------
void Resonance::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    DecayingParticle::addSymmetrizationIndex(c);
    MassShape_->addSymmetrizationIndex(c);
}

//-------------------------
void Resonance::clearSymmetrizationIndices()
{
    DecayingParticle::clearSymmetrizationIndices();
    MassShape_->clearSymmetrizationIndices();
}


}
