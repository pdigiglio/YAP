#include "Particle.h"

#include "logging.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, double m, std::string name) :
    std::enable_shared_from_this<Particle>(),
    QuantumNumbers_(q),
    Mass_(new RealParameter(m)),
    Name_(name)
{}

//-------------------------
void Particle::setMass(std::shared_ptr<RealParameter> m)
{
    if (!m)
        throw exceptions::Exception("mass is unset", "Particle::setMass");
    Mass_ = m;
}

//-------------------------
bool Particle::consistent() const
{
    bool C = true;

    if (Mass_->value() < 0.) {
        FLOG(ERROR) << "mass is negative";
        C &= false;
    }

    if (ParticleCombinations_.empty()) {
        FLOG(ERROR) << "ParticleCombinations_ is empty";
        C &= false;
    }

    return C;
}

//-------------------------
void Particle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    ParticleCombinations_.push_back(pc);
}

//-------------------------
void Particle::pruneParticleCombinations()
{
    if (!model())
        throw exceptions::Exception("Model not set", "DecayChannel::pruneParticleCombinations");

    // remove entries that don't trace back to the ISP
    for (auto it = ParticleCombinations_.begin(); it != ParticleCombinations_.end(); ) {
        // find the top-most parent
        auto pc = *it;
        while (pc->parent())
            pc = pc->parent();
        // check if it's not an ISP
        if (pc->indices().size() != model()->finalStateParticles().size())
            // erase
            it = ParticleCombinations_.erase(it);
        else
            it++;
    }

    if (ParticleCombinations_.empty())
        throw exceptions::Exception("ParticleCombinations empty after pruning", "Particle::pruneParticleCombinations");
}

//-------------------------
std::string to_string(const Particle& p)
{
    return p.name() + "(" + to_string(p.quantumNumbers()) + "), mass = " + std::to_string(p.mass()->value());
}


}
