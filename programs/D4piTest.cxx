#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "SpinUtilities.h"
#include "WignerD.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <assert.h>
#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP


int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);


    unsigned max2L(2 * 4);

    /// \todo Figure out clever way to find PDL file
    yap::ParticleFactory factory("evt.pdl");


    // initial state particle
    double radialSize = 1.;
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(421, radialSize);

    // final state particles
    yap::helicityStates piPlus = factory.createFinalStateParticle(211, {0, 2});
    yap::helicityStates piMinus = factory.createFinalStateParticle(-211, {1, 3});


    // rho rho
    yap::helicityStates rho = factory.createResonanceBreitWigner(113, radialSize);
    rho.addChannels(piPlus, piMinus, max2L);

    D->addChannels(rho, rho, max2L);

    // omega omega
    yap::helicityStates omega = factory.createResonanceBreitWigner(223, radialSize);
    omega.addChannels(piPlus, piMinus, max2L);

    D->addChannels(omega, omega, max2L);

    // rho omega
    D->addChannels(rho, omega, max2L);


    // a_1 channels
    yap::helicityStates sigma = factory.createResonanceBreitWigner(9000221, radialSize);
    sigma.addChannels(piPlus, piMinus, max2L);

    yap::helicityStates a_1 = factory.createResonanceBreitWigner(20213, radialSize);
    a_1.addChannels(sigma, piPlus, max2L);

    a_1.addChannels(rho, piPlus, max2L);

    D->addChannels(a_1, piMinus, max2L);


    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.createResonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);




    // consistency and optimizations
    assert(D->prepare());
    std::cout << "consistent! \n";

    // print stuff
    yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\nD symmetrizations: \n";
    for (std::shared_ptr<yap::ParticleCombination>& pc : D->particleCombinations())
        std::cout << std::string(*pc) << "\n";
    std::cout << "\n";

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices: \n";
    for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices: \n";
    for (auto& pc : D->helicityAngles().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->helicityAngles().symmetrizationIndex(pc) << "\n";

    D->printDecayChain();
    std::cout << "\n";

    D->printSpinAmplitudes();
    D->printDataAccessors();



    // test helicity angles
    TLorentzVector P(0.0, 0.0, 0.0, D->mass());
    Double_t masses[4] = { piPlus[0]->mass(), piMinus[0]->mass(), piPlus[0]->mass(), piMinus[0]->mass() };

    for (unsigned int iEvt = 0; iEvt < 1000; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, 4, masses);
        event.Generate();

        std::vector<TLorentzVector> momenta;
        for (unsigned i = 0; i < 4; ++i)
            momenta.push_back(*event.GetDecay(i));

        assert(D->addDataPoint(momenta));
    }

    // to test amplitude calculation, set all free amps to 1
    std::vector<yap::Amp> freeAmps = D->freeAmplitudes();
    for (yap::Amp& a : freeAmps)
        a = yap::Complex_1;
    D->setFreeAmplitudes(freeAmps);

    D->logLikelihood();

    /*for (yap::Amp& a : freeAmps)
        a *= 0.5;
    D->setFreeAmplitudes(freeAmps);

    D->logLikelihood();*/





    std::cout << "alright! \n";
}
