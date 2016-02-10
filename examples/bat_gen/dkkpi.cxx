// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "dkkpi.h"

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include "ParticleCombination.h"
#include "QuantumNumbers.h"
#include "Resonance.h"

#include <complex>

std::unique_ptr<yap::Model> dkkpi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{

    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // final state particles
    auto kPlus  = F.fsp(+321);
    auto kMinus = F.fsp(-321);
    auto piPlus = F.fsp(+211);

    auto M = std::make_unique<yap::Model>(std::move(SAC));

    M->setFinalState({kPlus, kMinus, piPlus});

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // phi
    // auto phi = std::make_shared<yap::Resonance>(F.quantumNumbers("phi"), 1010.e-3, "phi", radialSize, std::make_shared<yap::BreitWigner>());
    auto phi = std::make_shared<yap::Resonance>(yap::QuantumNumbers(2, 0), 1310.e-3, "phi", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(phi->massShape())->width()->setValue(20e-3);
    phi->addChannel({kPlus, kMinus});
    D->addChannel({phi, piPlus});

    std::cout << *phi << std::endl;

    /*
    // X_2
    auto X_2 = std::make_shared<yap::Resonance>(yap::QuantumNumbers(4, 0), 1.2, "X_2", radialSize, std::make_shared<yap::BreitWigner>());
    std::static_pointer_cast<yap::BreitWigner>(X_2->massShape())->width()->setValue(80e-3);
    X_2->addChannel({piPlus, kMinus});
    D->addChannel({X_2, kPlus});
    */

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = M->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    return M;
}
