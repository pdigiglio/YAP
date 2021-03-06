#include "Flatte.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DataPoint.h"
#include "DataPartition.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "Resonance.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
void Flatte::addChannel(std::shared_ptr<RealParameter> coupling, std::shared_ptr<RealParameter> mass)
{
    if (!coupling)
        throw exceptions::Exception("Coupling is unset", "Flatte::addChannel");
    if (!mass)
        throw exceptions::Exception("Mass is unset", "Flatte::addChannel");
    FlatteChannels_.push_back(FlatteChannel(coupling, mass));

    addParameter(coupling);
    addParameter(mass);
}

//-------------------------
void Flatte::addChannel(double coupling, double mass)
{
    addChannel(std::make_shared<RealParameter>(coupling), std::make_shared<RealParameter>(mass));
}

//-------------------------
void Flatte::calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const
{
    // precalculate
    std::vector<std::complex<double> > fc_c2;
    std::vector<std::complex<double> > fc_4m2c2;
    fc_c2.reserve(FlatteChannels_.size());
    fc_4m2c2.reserve(FlatteChannels_.size());
    for (const auto& fc : FlatteChannels_) {
        fc_4m2c2.push_back(Complex_1 * 4. * pow(fc.Mass->value() * fc.Coupling->value(), 2));
        fc_c2.push_back(Complex_1 * pow(fc.Coupling->value(), 2));
    }

    const double M2 = pow(mass()->value(), 2);

    for (auto& d : D) {

        const double m2 = model()->fourMomenta()->m2(d, pc);

        // calculate width term
        auto w = Complex_0;

        // sum of coupling * complex-breakup-momentum * 2 * i / m
        for (size_t i = 0; i < FlatteChannels_.size(); ++i)
            w += std::sqrt(fc_4m2c2[i] / m2 - fc_c2[i]);

        // T = 1 / (M^2 - m^2 - width-term)
        T()->setValue(1. / (M2 - m2 - w), d, si, D);
    }
}

//-------------------------
bool Flatte::consistent() const
{
    bool C = MassShapeWithNominalMass::consistent();

    for (const auto& fc : FlatteChannels_) {
        if (fc.Coupling->value() <= 0) {
            FLOG(ERROR) << "coupling constant <= 0";
            C &= false;
        }
        if (fc.Mass->value() <= 0) {
            FLOG(ERROR) << "mass <= 0";
            C &= false;
        }
    }

    return C;
}

}




