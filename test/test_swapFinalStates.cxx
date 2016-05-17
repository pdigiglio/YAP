#include <catch.hpp>
#include <catch_capprox.hpp>

#include <AmplitudePair.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <ZemachFormalism.h>

#include <cmath>

/**
 * Test that the amplitude remains the same after swapping the order of the final state particles
 */

yap::MassAxes populate_model(yap::Model& M, const yap::ParticleFactory& F, const std::vector<int>& FSP)
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
    D->addChannel({piK0, kPlus})->freeAmplitudes()[0]->setValue(0.5 * yap::Complex_1);

    auto piK1 = yap::Resonance::create(yap::QuantumNumbers(2, 0), 1.00, "piK1", 3., std::make_shared<yap::BreitWigner>(0.025));
    piK1->addChannel({piPlus, kMinus});
    D->addChannel({piK1, kPlus})->freeAmplitudes()[0]->setValue(1. * yap::Complex_1);

    auto piK2 = yap::Resonance::create(yap::QuantumNumbers(4, 0), 1.25, "piK2", 3., std::make_shared<yap::BreitWigner>(0.025));
    piK2->addChannel({piPlus, kMinus});
    D->addChannel({piK2, kPlus})->freeAmplitudes()[0]->setValue(30. * yap::Complex_1);

    return M.massAxes({{i_piPlus, i_kMinus}, {i_kMinus, i_kPlus}});
}

std::complex<double> calculate_model(yap::Model& M, const yap::MassAxes& A, std::vector<double> m2, yap::DataSet& data)
{
    // calculate four-momenta
    auto P = M.calculateFourMomenta(A, m2);

    // if failed, outside phase space
    if (P.empty())
        return std::numeric_limits<double>::quiet_NaN();

    // reset data set
    data = M.createDataSet();
    // add point
    data.add(P);

    // return amplitude
    return M.amplitude(data[0], data);
}

TEST_CASE( "swapFinalStates" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

//    auto F = yap::ParticleFactory((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    yap::ParticleFactory F("../../data/evt.pdl");

    // create models
    std::vector<std::unique_ptr<yap::Model> > Z;     // Zemach
    std::vector<yap::MassAxes> mZ; // always (pi+, K-), (K-, K+)
    std::vector<yap::DataSet> dZ;

    std::vector<std::unique_ptr<yap::Model> > H;     // Helicity
    std::vector<yap::MassAxes> mH; // always (pi+, K-), (K-, K+)
    std::vector<yap::DataSet> dH;

    std::vector<int> FSP = {F.pdgCode("K-"), F.pdgCode("pi+"), F.pdgCode("K+")};
    std::sort(FSP.begin(), FSP.end());
    do {

        // Zemach
        Z.emplace_back(new yap::Model(std::make_unique<yap::ZemachFormalism>()));
        mZ.push_back(populate_model(*Z.back(), F, FSP));
        dZ.push_back(Z.back()->createDataSet(1));

        // Helicity
        H.emplace_back(new yap::Model(std::make_unique<yap::HelicityFormalism>()));
        mH.push_back(populate_model(*H.back(), F, FSP));
        dH.push_back(H.back()->createDataSet(1));

    } while (std::next_permutation(FSP.begin(), FSP.end()));

    // get piK and KK mass ranges
    auto m2_piK_range = Z[0]->massRange(mZ[0][0]);
    auto m2_KK_range  = Z[0]->massRange(mZ[0][1]);

    const unsigned N = 30;
    // loop over phase space
    for (double m2_piK = m2_piK_range[0]; m2_piK <= m2_piK_range[1]; m2_piK += (m2_piK_range[1] - m2_piK_range[0]) / N) {
        for (double m2_KK = m2_KK_range[0]; m2_KK <= m2_KK_range[1]; m2_KK += (m2_KK_range[1] - m2_KK_range[0]) / N) {

            std::vector<std::complex<double> > amps_Z(Z.size(), 0.);
            std::vector<std::complex<double> > amps_H(H.size(), 0.);

            for (size_t i = 0; i < Z.size(); ++i) {
                FDEBUG("Calculate Zemach for FinalState combination " << i);
                amps_Z[i] = calculate_model(*Z[i], mZ[i], {m2_piK, m2_KK}, dZ[i]);

                FDEBUG("Calculate Helicity for FinalState combination " << i);
                amps_H[i] = calculate_model(*H[i], mH[i], {m2_piK, m2_KK}, dH[i]);
            }



            // print
            std::cout << m2_piK << ", " << m2_KK << " is " << ((std::isnan(real(amps_Z[0]))) ? "out" : "in") << " phase space" << std::endl;

            std::cout << "Zemach:                        Helicity:" << std::endl;
            //for (size_t i = 0; i < amps_Z.size(); ++i)
            for (size_t i = 0; i < 1; ++i) {
                double phaseDiff = arg(amps_Z[i]) - arg(amps_H[i]);

                std::cout << amps_Z[i] << " " << norm(amps_Z[i]) << "     " << amps_H[i] << " " << norm(amps_H[i])
                          << "      ratio Z/H = " <<  norm(amps_Z[i]) / norm(amps_H[i])
                          << "      rel. phase = " << phaseDiff / yap::rad_per_deg<double>() << "°"
                          << std::endl;
            }
            /*
            for (size_t i = 0; i < H.size(); ++i) {

                auto PCs = H[i]->helicityAngles()->particleCombinations();

                // print fsp 4 momenta
                for (auto& pc : PCs) {
                    std::cout << "FSP 4mom for " << to_string(*pc) << " = "
                            << to_string(H[i]->fourMomenta()->p(dH[i][0], pc)) << "\n";
                }

                std::cout << "piK = " << *mH[i][0];

                for (size_t j = 0; j < PCs.size(); ++j) {
                    std::cout << "\t :: " << *PCs[j] << ": "
                              << "("  << H[i]->helicityAngles()->phi(dH[i][0], PCs[j]) * yap::deg_per_rad<double>()
                              << ", " << H[i]->helicityAngles()->theta(dH[i][0], PCs[j]) * yap::deg_per_rad<double>()
                              << ")";
                }
                std::cout << std::endl;
            }*/

            // if nan, check that all are nan
            if (amps_Z[0] != amps_Z[0]) {
                for (size_t i = 1; i < amps_Z.size(); ++i)
                    REQUIRE ( amps_Z[i] != amps_Z[i] );
                for (size_t i = 0; i < amps_H.size(); ++i)
                    REQUIRE ( amps_H[i] != amps_H[i] );
            }
            // otherwise check that amplitudes are approximately the same
            else {
                // // check equality for Zemach
                for (size_t i = 1; i < amps_Z.size(); ++i)
                    REQUIRE ( amps_Z[i - 1] == Catch::Detail::CApprox( amps_Z[i] ) );

                // check equality for Helicity
                for (size_t i = 1; i < amps_H.size(); ++i)
                    REQUIRE ( amps_H[i - 1] == Catch::Detail::CApprox( amps_H[i] ) );

                // todo check if Zemach and Helicity have the same phase
                double phaseDiff = arg(amps_Z[0]) - arg(amps_H[0]);
                phaseDiff /= yap::rad_per_deg<double>();
                /*// 180 degree are ok????????????
                if (phaseDiff >= 180)
                    phaseDiff -= 180;
                if (phaseDiff <= -180)
                    phaseDiff += 180;*/
                //REQUIRE( phaseDiff == Approx(0) );

                /// \todo check if Zemach and Helicity have same magnitude
                //REQUIRE( norm(amps_Z[i]) == Approx(norm(amps_H[i])) );
            }

        }
    }

}
