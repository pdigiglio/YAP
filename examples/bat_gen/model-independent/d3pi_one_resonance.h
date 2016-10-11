/**
 *
 *    @file  d3pi_one_resonance.h
 *   @brief  
 *
 *    @date  09/30/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  __BAT__D3PI__RHO__ONLY__H
#define  __BAT__D3PI__RHO__ONLY__H

#include "Model.h"

#include "../bat_fit.h"

#include <memory>
#include <vector>

/// Creates a decay \f$D^+ \to \pi^+\pi^-\pi^+\f$ with only one resonance.
std::unique_ptr<yap::Model> d3pi_one_resonance(std::unique_ptr<yap::Model> M, const std::string& res_name);

/// Creates a decay with mass-bin resonances. For fitting.
std::unique_ptr<yap::Model> d3pi_one_resonance_binned(std::unique_ptr<yap::Model> M, const std::string& res_name);

bat_fit d3pi_one_resonance_fit(const std::string& name, std::unique_ptr<yap::Model> M, std::vector<std::vector<unsigned>> pcs = {});

#endif
