#include "ShowerIndex.h"

using namespace Herwig;


PersistentOStream& operator<<(PersistentOStream & os, const ShowerIndex & x) {
  os << x.id << x.interaction << x.timeFlag;
  return os;
}

PersistentIStream& operator>>(PersistentIStream & is, ShowerIndex & x) {
  int interactionInt = -1;
  int timeFlagInt = -1;
  is >> x.id >> interactionInt >> timeFlagInt;
  x.interaction = ShowerIndex::InteractionType(interactionInt);
  x.timeFlag    = ShowerIndex::TimeOrderType(timeFlagInt);
  return is;
}

bool ShowerIndex::operator< (const ShowerIndex & rhs) const {
  bool order = true;
  if ( this->id > rhs.id ) {
    order = false;
  } else if ( this->id == rhs.id ) {
    if ( this->interaction > rhs.interaction ) {
      order = false;
    } else if ( this->interaction == rhs.interaction ) { 
      if ( this->timeFlag >= rhs.timeFlag ) {
	order = false;
      }
    }
  }
  return order;
}


