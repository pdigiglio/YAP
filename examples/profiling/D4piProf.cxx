/**
 *
 *    @file  D4piProf.cxx
 *   @brief  Simple example to profile \f$D\to\pi^+\pi^+\pi^-\pi^-\f$
 *
 *    @date  06/08/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */
#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "Flatte.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "FreeAmplitude.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "Model.h"
#include "Parameter.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PoleMass.h"
#include "Resonance.h"
#include "SpinAmplitudeCache.h"
#include "WignerD.h"
#include "logging.h"
#include "make_unique.h"

#include <iostream>
int main (int argc, char *argv[]) {

	yap::Model M(std::make_unique<yap::HelicityFormalism>());

	yap::ParticleFactory f((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

	// Set initial state particle
	const double radialSize = 1.; // [1/GeV]
	auto D = f.decayingParticle(421, radialSize);

	// Set final state particles
	auto piPlus  = f.fsp(211);
	auto piMinus = f.fsp(-211);

	// Set final state
	M.setFinalState({piPlus, piMinus, piPlus, piMinus});

    // sigma
    auto sigma = f.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel({piPlus, piMinus});

    // rho
    auto rho = f.resonance(113, radialSize, std::make_shared<yap::PoleMass>());
    rho->addChannel({piPlus, piMinus});

    // omega
    auto omega = f.resonance(223, radialSize, std::make_shared<yap::BreitWigner>());
    omega->addChannel({piPlus, piMinus});

    // a_1
	auto a_1_flatte = std::make_shared<yap::Flatte>();
	a_1_flatte->addChannel(.5, .5);
	a_1_flatte->addChannel(.5, .5);
    auto a_1 = f.resonance(20213, radialSize, a_1_flatte);
    a_1->addChannel({sigma, piPlus});
    a_1->addChannel({rho,   piPlus});

    // D's channels
    D->addChannel({rho, rho});
    D->addChannel({omega, omega});
    D->addChannel({rho, omega});
    D->addChannel({a_1, piMinus});


	if (!M.consistent()) {
		std::cerr << "Inconsistent. Exiting." << std::endl;
		return 1;
	}

	return 0;
}
