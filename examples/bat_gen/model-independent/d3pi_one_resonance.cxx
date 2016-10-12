/**
 *
 *    @file  d3pi_one_resonance.cxx
 *   @brief  
 *
 *    @date  09/30/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "./d3pi_one_resonance.h"
#include "./MassBin.h"

#include "../bat_fit.h"
#include "../ConstantPrior.h"
#include "../fit_fitFraction.h"

#include <BreitWigner.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <Model.h>
#include <PDL.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>
#include <container_utils.h>
#include <make_unique.h>

//#include <BAT/BCGaussianPrior.h>
//#include <BAT/BCSplitGaussianPrior.h>

#include <complex>
#include <iterator>
#include <memory>

std::vector<double> make_partition(double min, const double max, const unsigned int bins)
{
    std::vector<double> partition;
    partition.reserve(bins + 1);

    const double step = (max - min) / bins;
    for (unsigned int i = 0; i < bins; ++i) {
        std::cout << "min = " << min << std::endl;
        partition.push_back(min);
        min += step;
    }
    // insert last element
    partition.push_back(max);

    return partition;
}

std::unique_ptr<yap::Model> d3pi_binned(std::unique_ptr<yap::Model> M, const std::string& res_name)
{
    using namespace std;
    using namespace yap;

    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = F.fsp(211);
    auto piMinus = F.fsp(-211);

    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // number of bins
    constexpr unsigned bins = 20;
    // partition
    auto p = make_partition(std::abs(piPlus->mass() + piMinus->mass()),
                            std::abs(F.fsp(F.pdgCode("D+"))->mass() - piPlus->mass()), bins);
    // vector of resonance pointers
    std::vector<shared_ptr<Resonance>> res;
    res.reserve(bins);
    for (int i = 0; i < static_cast<int>(p.size()) - 1; ++i) {
        auto r = Resonance::create(
                res_name + "_" + to_string(i),
                QuantumNumbers(0, 0),
                radialSize,
                make_shared<MassBin>(p[i], p[i+1]));
        std::cout << "bin[" << i << "] = " << p[i] << " to " << p[i+1] << std::endl;
        r->addChannel(piPlus, piMinus);
        D->addChannel(r, piPlus);
//        res.emplace_back(move(r));
    }

    M->addInitialStateParticle(D);

    
    auto odd_numbers = [](int i) { return i % 2; };
    int i = 0;
    std::cout << ">>>>>>>>>>> bins = " << bins
              << " free_amplitudes = " << free_amplitudes(*M, from(*D)).size()
              << std::endl;

    auto delta_function = [](int i, int j) { return static_cast<double>(i == j); };
    for (const auto& fa : free_amplitudes(*M, from(*D))) {
        *fa = polar( i & 1 ? 1. : 10.,  0.);
        std::cout << "i = " << i << std::endl;
        std::cout << to_string(*fa) << std::endl;
        ++i;
    }

    return M;
}

std::unique_ptr<yap::Model> d3pi_one_resonance(std::unique_ptr<yap::Model> M, const std::string& res_name)
{
    using namespace std;
    using namespace yap;

    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = F.fsp(211);
    auto piMinus = F.fsp(-211);

    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    auto r = F.resonance(F.pdgCode(res_name), radialSize, make_shared<BreitWigner>());
    r->addChannel(piPlus, piMinus);
    D->addChannel(r, piPlus);

    M->addInitialStateParticle(D);

    *free_amplitude(*M, to(r)) = polar(1., 0.);

    return M;
}

bat_fit d3pi_binned_fit(const std::string& name,
                               std::unique_ptr<yap::Model> M,
                               std::vector<std::vector<unsigned>> pcs)
{
    using namespace std;
    using namespace yap;

    bat_fit m(name, d3pi_binned(move(M), "f_0"), pcs);

    m.model()->lock();

    for (const auto& fa : free_amplitudes(*m.model(), is_not_fixed())) {
        cout << "Fit: fa = " << to_string(*fa) << std::endl;
        m.setPriors(fa, new ConstantPrior(0, 11), new ConstantPrior(-180, 180));
        m.setRealImagRanges(fa, -11, 11, -11, 11);
        m.setAbsArgRanges(fa, 0, 11, -180, 180);
    }

    // Fix an amplitude
    // Use the scope to have local variables
    { 
        const auto fas = free_amplitudes(*m.model(), is_not_fixed());
        const auto& fa = next(fas.cbegin(), 1);
        if (fa == fas.cend())
            throw 1;

        std::cout << "Fit: free amplitudes = " << fas.size() << std::endl;
        m.fix(*fa, real(1.), imag(0.));
    }

    return m;
}
