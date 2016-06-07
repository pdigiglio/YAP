// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKKPI__H
#define __BAT__DKKPI__H

#include "BreitWigner.h"
#include "Constants.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include <Model.h>
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "QuantumNumbers.h"
#include "Resonance.h"
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

inline std::unique_ptr<yap::Model> dkkpi(std::unique_ptr<yap::SpinAmplitudeCache> SAC)
{
    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

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

    auto KK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 1.1, "KK0", radialSize, std::make_shared<yap::BreitWigner>(0.025));
    KK0->addChannel({kPlus, kMinus});
    D->addChannel({KK0, piPlus})->freeAmplitudes().begin()->get()->setValue(1. * yap::Complex_1);

    auto KK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.35, "KK1", radialSize, std::make_shared<yap::BreitWigner>(0.025));
    KK1->addChannel({kPlus, kMinus});
    D->addChannel({KK1, piPlus})->freeAmplitudes().begin()->get()->setValue(std::polar<double>(2, 0 * yap::rad_per_deg<double>()));

    auto KK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.6, "KK2", radialSize, std::make_shared<yap::BreitWigner>(0.025));
    KK2->addChannel({kPlus, kMinus});
    D->addChannel({KK2, piPlus})->freeAmplitudes().begin()->get()->setValue(30. * yap::Complex_1);

    auto piK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.75, "piK0", radialSize, std::make_shared<yap::BreitWigner>(0.025));
    piK0->addChannel({piPlus, kMinus});
    D->addChannel({piK0, kPlus})->freeAmplitudes().begin()->get()->setValue(std::polar<double>(1, 180 * yap::rad_per_deg<double>()));

    auto piK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.00, "piK1", radialSize, std::make_shared<yap::BreitWigner>(0.025));
    piK1->addChannel({piPlus, kMinus});
    D->addChannel({piK1, kPlus})->freeAmplitudes().begin()->get()->setValue(1. * yap::Complex_1);

    // auto piK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.25, "piK2", radialSize, std::make_shared<yap::BreitWigner>(1.25, 0.025));
    // piK2->addChannel({piPlus, kMinus});
    // D->addChannel({piK2, kPlus})->freeAmplitudes().begin()->get()->setValue(1. * yap::Complex_1);

    return M;
}

#endif
