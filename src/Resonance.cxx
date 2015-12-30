#include "Resonance.h"

#include "DecayChannel.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombinationCache.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::unique_ptr<MassShape> massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(std::move(massShape))
{
    if (!MassShape_)
        throw exceptions::MissingMassShape();

    MassShape_->setResonance(this);
}

//-------------------------
bool Resonance::consistent() const
{
    bool C = DecayingParticle::consistent();

    if (!MassShape_->consistent()) {
        FLOG(ERROR) << "MassShape not consistent";
        C &= false;
    }
    if (MassShape_->resonance() != this) {
        FLOG(ERROR) << "MassShape's resonance is not this";
        C &= false;
    }

    return C;
}

//-------------------------
void Resonance::addChannel(std::unique_ptr<DecayChannel> c)
{
    if (!initialStateParticle())
        throw exceptions::InitialStateParticleUnset();

    for (auto& pc : c->particleCombinations())
        MassShape_->addSymmetrizationIndex(initialStateParticle()->particleCombinationCache[pc]);

    DecayingParticle::addChannel(std::move(c));
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
