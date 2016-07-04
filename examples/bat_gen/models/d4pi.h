// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <QuantumNumbers.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

inline std::unique_ptr<yap::Model> d4pi()
{
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    M->setFinalState({piPlus, piMinus, piPlus, piMinus});

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);
    
    //
    // resonant particles
    //
    
    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel({piPlus, piMinus});

    // omega
    //auto omega = F.resonance(F.pdgCode("omega"), radialSize, std::make_shared<yap::BreitWigner>());
    //omega->addChannel({piPlus, piMinus});
    
    // sigma / f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel({piPlus, piMinus});

    // a_1
    auto a_1 = F.resonance(F.pdgCode("a_1+"), radialSize, std::make_shared<yap::BreitWigner>());
    
    for (auto& freeAmp : a_1->addChannel({rho, piPlus})->freeAmplitudes()) {
        LOG(INFO) << to_string(*freeAmp);
        if (freeAmp->spinAmplitude()->L() == 0)
            freeAmp->setVariableStatus(yap::VariableStatus::fixed); // S-wave, fixed
        else if (freeAmp->spinAmplitude()->L() == 1)
            freeAmp->setValue(0.); // P-wave, 0
        else if (freeAmp->spinAmplitude()->L() == 2)
            freeAmp->setValue(std::polar(0.241, yap::rad(82.))); // D-wave
        else
            LOG(ERROR) << "wrong spin";
    }
    
    a_1->addChannel({sigma, piPlus})->freeAmplitudes().begin()->get()->setValue(std::polar(0.439, yap::rad(193.)));
    
    // f_0(980) (as Flatte)
    auto piZero = F.fsp(111);
    auto Kshort = F.fsp(310);
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.20, piZero->mass()->value());
    f_0_980_flatte->addChannel(0.50, Kshort->mass()->value());
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addChannel({piPlus, piMinus});
       
    // f_2(1270)
    auto f_2 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<yap::BreitWigner>());
    f_2->addChannel({piPlus, piMinus}); 
    
    // pi+ pi- flat
    auto pipiFlat = F.decayingParticle(F.pdgCode("f_0"), radialSize); // just need any spin0 particle
    pipiFlat->addChannel({piPlus, piMinus});   
    
    //
    // D0 channels
    //
    
    D->addChannel({a_1, piMinus});
    
    
    // rho rho
    // transform into angular momentum basis
    // A_S = sqrt(6/5) * A_parallel - sqrt(3/5) * A_0
    // A_P = A_perp
    // A_D = sqrt(3/5) * A_parallel + sqrt(6/5) * A_0
    auto A_parallel = std::polar(0.157, yap::rad(120.));
    auto A_perp     = std::polar(0.384, yap::rad(163.));
    auto A_0        = std::polar(0.624, yap::rad(357.));
    
    for (auto& freeAmp : D->addChannel({rho, rho})->freeAmplitudes()) {
        LOG(INFO) << to_string(*freeAmp);
        if (freeAmp->spinAmplitude()->L() == 0)
            freeAmp->setValue(sqrt(6./5.) * A_parallel - sqrt(3./5.) * A_0); // S-wave
        else if (freeAmp->spinAmplitude()->L() == 1)
            freeAmp->setValue(A_perp); // P-wave
        else if (freeAmp->spinAmplitude()->L() == 2)
            freeAmp->setValue(sqrt(3/5) * A_parallel + sqrt(6/5) * A_0); // D-wave
        else
            LOG(ERROR) << "wrong spin";
    }
    
    // R pi pi
    D->addChannel({f_0_980, piPlus, piMinus})->freeAmplitudes().begin()->get()->setValue(std::polar(0.233, yap::rad(261.)));
    D->addChannel({f_2,     pipiFlat       })->freeAmplitudes().begin()->get()->setValue(std::polar(0.338, yap::rad(317.)));
    D->addChannel({sigma,   piPlus, piMinus})->freeAmplitudes().begin()->get()->setValue(std::polar(0.432, yap::rad(254.)));

    M->addInitialStateParticle(D);

    return M;
}

#endif
