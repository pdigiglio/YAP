/// \file

#ifndef  yap_MassWindow_h
#define  yap_MassWindow_h

#include "fwd/DataPartition.h"
#include "fwd/ParticleCombination.h"
#include "fwd/Parameter.h"

#include "MassShapeWithNominalMass.h"

#include <memory>

/// Class to model a mass bin \f$m_{\text{min}} \le m < m_{\text{max}}\f$.
class MassBin : public yap::MassShapeWithNominalMass
{
public:
    /// \brief Constructor.
    /// \param lower_edge The lower value of the mass.
    /// \param upper_edge The upper value of the mass.
    explicit MassBin(double lower_edge, double upper_edge);

    double lowerEdge() const noexcept
    { return LowerEdge_; }

    double upperEdge() const noexcept
    { return UpperEdge_; }

private:
    /// Lower edge of the bin.
    const double LowerEdge_;
    /// Upper edge of the bin.
    const double UpperEdge_;

protected:

    virtual void calculateT(yap::DataPartition& D,
                            const std::shared_ptr<const yap::ParticleCombination>& pc,
                            unsigned si) const override;

    /// Dumb parameter to have the integration triggered.
//    std::shared_ptr<yap::RealParameter> DumbParam_;

};

#endif
