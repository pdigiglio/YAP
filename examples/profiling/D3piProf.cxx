/**
 *
 *    @file  profiling/D3piProf.cxx 
 *   @brief  
 *
 *    @date  06/07/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */
#include "logging.h"
#include "BreitWigner.h"
//#include "DataSet.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
//#include "HelicityAngles.h"
//#include "HelicityFormalism.h"
//#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
//#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "ZemachFormalism.h"

#include <iostream>

int main (int argc, char *argv[]) {
	yap::Model M(std::make_unique<yap::ZemachFormalism>());

	const char* basePath = std::getenv("YAPDIR");
	yap::ParticleFactory f((basePath ? ( static_cast<std::string>(basePath) + "/data/" ) : ("")) + "evt.pdl"); 

	// Set initial state particle
	const double radialSize = 3.; // [1/GeV]
	auto D = f.decayingParticle(f.pdgCode("D+"), radialSize);

	// Set final state particles
	auto piPlus  = f.fsp(221);
	auto piMinus = f.fsp(-221);

	// Set final state
	M.setFinalState({piPlus, piMinus, piPlus});

	// Rho
	auto rho = yap::Resonance::create(f.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(rho->massShape())->width()->setValue(0.149);
	rho->addChannel({piPlus, piMinus});

	// f_2(1270)
	auto f_2 = yap::Resonance::create(f.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(f_2->massShape())->width()->setValue(0.185);
	f_2->addChannel({piPlus, piMinus});

	// f_0(980)
	auto f_0_980 = yap::Resonance::create(f.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(f_0_980->massShape())->width()->setValue(0.329);
	f_0_980->addChannel({piPlus, piMinus});

	// f_0(1370)
	auto f_0_1370 = yap::Resonance::create(f.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->width()->setValue(0.250);
	f_0_1370->addChannel({piPlus, piMinus});

	// f_0(1500)
	auto f_0_1500 = yap::Resonance::create(f.quantumNumbers("f_0"), 1.507, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(f_0_1500->massShape())->width()->setValue(0.109);
	f_0_1500->addChannel({piPlus, piMinus});

	// sigma a.k.a. f_0(500)
	auto sigma = yap::Resonance::create(f.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
	std::static_pointer_cast<yap::BreitWigner>(sigma->massShape())->width()->setValue(0.800);
	sigma->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});
    D->addChannel({f_2,      piPlus});
    D->addChannel({f_0_980,  piPlus});
    D->addChannel({f_0_1370, piPlus});
    D->addChannel({f_0_1500, piPlus});
    D->addChannel({sigma,    piPlus});

	return 0;
}
