#include "ShowerIndex.h"

using namespace Herwig;


PersistentOStream & operator<< (PersistentOStream & os, const ShowerIndex & x) {
  os << x.id << x.interaction << x.timeFlag;
  return os;
}

PersistentIStream & operator>> (PersistentIStream & is, ShowerIndex & x) {
  int interactionInt = -1;
  int timeFlagInt = -1;
  is >> x.id >> interactionInt >> timeFlagInt;
  x.interaction = ShowerIndex::int2Interaction( interactionInt );
  x.timeFlag    = ShowerIndex::int2TimeOrder( timeFlagInt );
  return is;
}


ShowerIndex::InteractionType ShowerIndex::int2Interaction(const int position) {
  ShowerIndex::InteractionType theInteraction;
  switch ( position ) {
  case 0:  theInteraction = ShowerIndex::QCD; break;
  case 1:  theInteraction = ShowerIndex::QED; break;
  case 2:  theInteraction = ShowerIndex::EWK; break;
  default: theInteraction = ShowerIndex::UNDEFINED;
  }
  return theInteraction;
}

ShowerIndex::TimeOrderType ShowerIndex::int2TimeOrder(const int position) {
  ShowerIndex::TimeOrderType theTimeFlag;
  switch ( position ) {
  case 0:  theTimeFlag = ShowerIndex::IS; break;
  case 1:  theTimeFlag = ShowerIndex::FS; break;
  default: theTimeFlag = ShowerIndex::UNINITIALIZED;
  }
  return theTimeFlag;
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


