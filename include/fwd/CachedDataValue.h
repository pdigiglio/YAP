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

#ifndef yap_CachedDataValueFwd_h
#define yap_CachedDataValueFwd_h

#include <memory>
#include <set>
#include <vector>

namespace yap {

class CachedDataValue;
struct DaughterCachedDataValue;
class RealCachedDataValue;
class ComplexCachedDataValue;
class FourVectorCachedDataValue;

// Cannot be forward declared since it is nested inside
// CachedDataValue
/* class CachedDataValue::Status; */

/// \typedef CachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueSet = std::set<std::shared_ptr<CachedDataValue> >;

/// \typedef CachedDataValueVector
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueVector = std::vector<std::shared_ptr<CachedDataValue> >;

/// \typedef DaughterCachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using DaughterCachedDataValueSet = std::set<DaughterCachedDataValue>;

}

#endif
