#include <catch.hpp>
#include <catch_capprox.hpp>

#include <BreitWigner.h>
#include <DataPartition.h>
#include <DataPoint.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>
#include <WignerD.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <make_unique.h>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <cmath>

void Populate_model(yap::Model& M, const yap::ParticleFactory& F, const std::vector<int>& FSP)
{
    // create and set final-state particles
    M.setFinalState({F.fsp(FSP[0]), F.fsp(FSP[1]), F.fsp(FSP[2])});

    // find FSP's
    unsigned i_piPlus = FSP.size();
    unsigned i_kPlus = FSP.size();
    unsigned i_kMinus = FSP.size();
    for (size_t i = 0; i < FSP.size(); ++i)
        if (FSP[i] == F.pdgCode("pi+"))
            i_piPlus = i;
        else if (FSP[i] == F.pdgCode("K+"))
            i_kPlus  = i;
        else if (FSP[i] == F.pdgCode("K-"))
            i_kMinus = i;
    auto piPlus = M.finalStateParticles()[i_piPlus];
    auto kPlus = M.finalStateParticles()[i_kPlus];
    auto kMinus = M.finalStateParticles()[i_kMinus];

    // create ISP
    auto D = F.decayingParticle(F.pdgCode("D+"), 3.);

    // create resonances
    auto piK0 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.75, "piK0", 3., std::make_shared<yap::BreitWigner>(0.025));
    piK0->addChannel({piPlus, kMinus});
    D->addChannel({piK0, kPlus})->freeAmplitudes().begin()->get()->setValue(0.5 * yap::Complex_1);

    auto piK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.00, "piK1", 3., std::make_shared<yap::BreitWigner>(0.025));
    piK1->addChannel({piPlus, kMinus});
    D->addChannel({piK1, kPlus})->freeAmplitudes().begin()->get()->setValue(1. * yap::Complex_1);

    auto piK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.25, "piK2", 3., std::make_shared<yap::BreitWigner>(0.025));
    piK2->addChannel({piPlus, kMinus});
    D->addChannel({piK2, kPlus})->freeAmplitudes().begin()->get()->setValue(30. * yap::Complex_1);

    M.massAxes({{i_piPlus, i_kMinus}, {i_kMinus, i_kPlus}});
}

TEST_CASE( "Integration" )
{
    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    std::vector<int> FSP = {F.pdgCode("K-"), F.pdgCode("pi+"), F.pdgCode("K+")};

    yap::Model M(std::make_unique<yap::ZemachFormalism>());
    Populate_model(M, F, FSP);
    auto D = M.initialStateParticle();

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[4];
    for (size_t i = 0; i < FSP.size(); ++i) {
        masses[i] = M.finalStateParticles()[i]->mass()->value();
    }

    LOG(INFO) << "create dataPoints";

    // create data set
    unsigned nPoints = 1e5;
    yap::DataSet data = M.createDataSet(nPoints);
    yap::DataSet dataTest = M.createDataSet(nPoints);

    for (unsigned int iEvt = 0; iEvt < nPoints; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, FSP.size(), masses);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        for (unsigned i = 0; i < FSP.size(); ++i) {
            TLorentzVector p = *event.GetDecay(i);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));

            DEBUG(yap::to_string(momenta.back()));
        }

        data[iEvt].setFinalStateMomenta(momenta);
        dataTest[iEvt].setFinalStateMomenta(momenta);
    }

	M.calculate(data);

	unsigned nChains = 2;
	auto parts = yap::DataPartitionWeave::create(data, nChains);

	for (unsigned i = 0; i < nChains; ++ i ) {
		double integralSum = 0;
		for (auto& dtv: M.initialStateParticle()->decayTrees()) {
			for (auto& dt: dtv.second)
				integralSum += yap::integral(*dt, *parts[i]);
		}

		LOG(INFO) << integralSum;
	}

//	for (auto& fa: yap::freeAmplitudes(M.initialStateParticle()->decayTrees())) {
//		LOG(INFO) << fa->value();
//	}
	
}
