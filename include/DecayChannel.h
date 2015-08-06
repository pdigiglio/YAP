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

#ifndef yap_DecayChannel_h
#define yap_DecayChannel_h

#include "DataAccessor.h"

#include "BlattWeisskopf.h"
#include "Particle.h"
#include "SpinAmplitude.h"
#include <tuple>

namespace yap {

typedef std::array<Particle*, 2> Daughters;

class DecayChannel : public DataAccessor {
public:
  DecayChannel(Particle* daughterA, Particle* daughterB, unsigned int L, SpinAmplitude& spinAmplitude);
  ~DecayChannel() {;}

  virtual Amp amplitude(DataPoint& d);
  virtual bool checkConsistency() const;

  const Daughters& getDaughters() const {return Daughters_;}
  const Particle* getDaughter(unsigned int i) {return Daughters_.at(i);}

private:
  Daughters Daughters_;
  unsigned int L_; /// relative angular momentum between daughters
  BlattWeisskopf BlattWeisskopf_;
  SpinAmplitude& SpinAmplitude_; /// SpinAmplitude can be shared between several DecayChannels
  Amp FreeAmplitude_;
};

}

#endif
