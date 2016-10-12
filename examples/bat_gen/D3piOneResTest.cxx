#include "./model-independent/MassBin.h"
#include "./model-independent/d3pi_one_resonance.h"
#include "tools.h"

#include "BreitWigner.h"
#include "container_utils.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "DecayTreeVectorIntegral.h"
#include "Filters.h"
#include "FinalStateParticle.h"
#include "Flatte.h"
#include "FourMomenta.h"
#include "FreeAmplitude.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "ImportanceSampler.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "MassRange.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PoleMass.h"
#include "PHSP.h"
#include "RelativisticBreitWigner.h"
#include "Resonance.h"
#include "ZemachFormalism.h"

#include <memory>
#include <random>
#include <vector>

int main (int argc, const char* argv[]) {
    yap::plainLogs(el::Level::Debug);
    auto M = d3pi_binned(yap_model<yap::ZemachFormalism>(), "f_0");

    constexpr double D_mass = 1.8696200e+00;
    // get default Dalitz axes
    auto A   = M->massAxes();
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M->finalStateParticles()));

    // generate points randomly in phase space of model
    std::mt19937 g(0);

    // create data set
    auto data = M->createDataSet();

    LOG(INFO) << "DataSet created" ;
    // generate 10,000 phase-space-distributed data points
    std::generate_n(std::back_inserter(data), 100,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*M), D_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO) << data.size() << " data points of " << data[0].bytes() << " bytes each = " << (data.size() * data[0].bytes()) * 1.e-6 << " MB";

    M->calculate(data);

    std::cout << "sum of log intensity = " << sum_of_log_intensity(*M, data, 0) << std::endl; 

    yap::ModelIntegral MI(*M);
    yap::ImportanceSampler::calculate(MI, data);

    for (const auto& b2_dtvi : MI.integrals()) {

        LOG(INFO) << "\n" << to_string(b2_dtvi.second.decayTrees());

        auto A_DT = amplitude(b2_dtvi.second.decayTrees(), data[0]);
        LOG(INFO) << "A_DT = " << A_DT;
        LOG(INFO) << "|A_DT|^2 = " << norm(A_DT);

        LOG(INFO) << "integral = " << to_string(integral(b2_dtvi.second));
        auto ff = fit_fractions(b2_dtvi.second);
        for (size_t i = 0; i < ff.size(); ++i)
            LOG(INFO) << "fit fraction " << to_string(100. * ff[i]) << "% for " << to_string(*b2_dtvi.second.decayTrees()[i]->freeAmplitude());
        LOG(INFO) << "sum of fit fractions = " << to_string(std::accumulate(ff.begin(), ff.end(), yap::RealIntegralElement()));

        LOG(INFO) << "cached integral components:";
        auto I_cached = cached_integrals(b2_dtvi.second);
        for (const auto& row : I_cached)
            LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                         [](std::string & s, const yap::ComplexIntegralElement & c)
        { return s += "\t" + to_string(c);}).erase(0, 1);

        LOG(INFO) << "integral components:";
        auto I = integrals(b2_dtvi.second);
        for (const auto& row : I)
            LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                         [](std::string & s, const yap::ComplexIntegralElement & c)
        { return s += "\t" + to_string(c);}).erase(0, 1);
    }

    LOG(INFO) << "alright!";

    return 0;
}
