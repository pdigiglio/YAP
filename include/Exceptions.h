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

#ifndef yap_Exceptions_h
#define yap_Exceptions_h

#include <stdexcept>

namespace yap {

namespace exceptions {

/// \class AngularMomentumNotConserved
/// \defgroup Exceptions YAP exceptions
class AngularMomentumNotConserved : public std::exception {};

/// \class MissingSpinAmplitude
/// \ingroup Exceptions
class MissingSpinAmplitude : public std::exception {};

/// \class MissingMassShape
/// \ingroup Exceptions
class MissingMassShape : public std::exception {};

/// \class MissingDecayChannel
/// \ingroup Exceptions
class MissingDecayChannel : public std::exception {};

/// \class IncorrectDaughters
/// \ingroup Exceptions
class IncorrectDaughters : public std::exception {};

/// \class NoDaughters
/// \ingroup Exceptions
class NoDaughters : public IncorrectDaughters {};

/// \class OnlyOneDaughter
/// \ingroup Exceptions
class OnlyOneDaughter : public IncorrectDaughters {};

/// \class FewerThanTwoDaughters
/// \ingroup Exceptions
class FewerThanTwoDaughters : public IncorrectDaughters {};

/// \class MoreThanTwoDaughters
/// \ingroup Exceptions
class MoreThanTwoDaughters : public IncorrectDaughters {};

/// \class EmptyDaughter
/// \ingroup Exceptions
class EmptyDaughter : public IncorrectDaughters {};

/// \class InvalidDaughter
/// \ingroup Exceptions
class InvalidDaughter : public std::exception {};

/// \class InitialStateParticleMismatch
/// \ingroup Exceptions
class InitialStateParticleMismatch : public std::exception {};

/// \class InitialStateParticleUnset
/// \ingroup Exceptions
class InitialStateParticleUnset : public std::exception {};

/// \class DecayChannelEmpty
/// \ingroup Exceptions
class DecayChannelEmpty : public std::exception {};

/// \class DecayChannelAlreadySet
/// \ingroup Exceptions
class DecayChannelAlreadySet : public std::exception {};

/// \class ParticleCombinationsEmpty
/// \ingroup Exceptions
class ParticleCombinationsEmpty : public std::exception {};

/// \class ParticleCombinationNotFound
/// \ingroup Exceptions
class ParticleCombinationNotFound : public std::exception {};

/// \class InitialStateParticleAlreadyPrepared
/// \ingroup Exceptions
class InitialStateParticleAlreadyPrepared : public std::exception {};

/// \class Inconsistent
/// \ingroup Exceptions
class Inconsistent : public std::exception {};

/// \class InconsistentSpinProjection
/// \ingroup Exceptions
class InconsistentSpinProjection : public Inconsistent {};

/// \class InconsistentParticleCombinationCache
/// \ingroup Exceptions
class InconsistentParticleCombinationCache : public Inconsistent {};

}

}

#endif
