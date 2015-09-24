/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_DataPoint_h
#define yap_DataPoint_h

#include "Amp.h"

#include <TLorentzVector.h>

#include <memory>
#include <set>
#include <vector>

namespace yap {

class DataAccessor;
class FourMomenta;
class HelicityAngles;
class ParticleCombination;

/// \class DataPoint
/// \brief Class for holding data and cached values per data point for fast calculation
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Data Data-related classes

class DataPoint
{
public:

    /// \name Constructors
    /// @{

    /// 4-momenta constructor
    DataPoint(const std::vector<TLorentzVector>& P);

    // /// Invariant mass constructor
    // DataPoint(const std::vector<double>& S);

    /// @}

    bool setFourMomenta(const std::vector<TLorentzVector>& fourMomenta);

    void printDataSize();

    /// \name Data accessor friends
    /// @{

    friend class FourMomenta;
    friend class HelicityAngles;
    friend class DataAccessor;
    friend class DataSet;

    /// reserve space in vectors
    void allocateStorage(const FourMomenta& fourMom, const HelicityAngles& helAngles, const std::set<DataAccessor*> dataAccessors);

    /// @}

protected:

    /// Vector of 4-momenta of particles in event
    std::vector<TLorentzVector> FourMomenta_;

    /// Helicity angles phi and theta
    /// first index is for the symmeterization state (as known by the InitialStateParticle)
    /// second index is for [phi, theta]
    std::vector<std::vector<double> > HelicityAngles_;

    /// Data storage for all DataAccessors

    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<std::vector<double> > > Data_;

    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    std::vector<std::vector<Amp> > CachedAmplitudes_;

};

}

#endif
