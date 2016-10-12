#include "./MassBin.h"

#include "DataPartition.h"
#include "FourMomenta.h"
#include "Model.h"
#include "Parameter.h"

#include <complex>
#include <memory>

//-------------------------
MassBin::MassBin(double lower_edge, double upper_edge) :
    yap::MassShapeWithNominalMass(),
    LowerEdge_(lower_edge),
    UpperEdge_(upper_edge)
//    DumbParam_(std::make_shared<yap::RealParameter>(0.))
{
//    addParameter(DumbParam_);
}

//-------------------------
void MassBin::calculateT(yap::DataPartition& D,
                         const std::shared_ptr<const yap::ParticleCombination>& pc,
                         unsigned si) const
{
    auto in_range = [=](double m) { return LowerEdge_ <= m && m < UpperEdge_; };
    for (auto& d : D) {
        auto   mass        = model()->fourMomenta()->m(d, pc);
        double is_in_range = in_range(mass);

//        std::cout << mass << " is between " << LowerEdge_
//                  << " and " << UpperEdge_
//                  << " = "   << is_in_range << std::endl;

        T()->setValue(std::complex<double>(1. * is_in_range, 0), d, si, D);
    }
}
