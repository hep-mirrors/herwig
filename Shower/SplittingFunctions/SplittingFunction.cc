// -*- C++ -*-
//
// SplittingFunction.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingFunction class.
//

#include "SplittingFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

AbstractClassDescription<SplittingFunction> SplittingFunction::initSplittingFunction;
// Definition of the static class description member.

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("The SplittingFunction class is the based class for 1->2 splitting functions"
     " in Herwig++");

  static Switch<SplittingFunction,ColourStructure> interfaceColourStructure
    ("ColourStructure",
     "The colour structure for the splitting function",
     &SplittingFunction::_colourStructure, Undefined, false, false);
  static SwitchOption interfaceColourStructureTripletTripletOctet
    (interfaceColourStructure,
     "TripletTripletOctet",
     "3 -> 3 8",
     TripletTripletOctet);
  static SwitchOption interfaceColourStructureOctetOctetOctet
    (interfaceColourStructure,
     "OctetOctetOctet",
     "8 -> 8 8",
     OctetOctetOctet);
  static SwitchOption interfaceColourStructureOctetTripletTriplet
    (interfaceColourStructure,
     "OctetTripletTriplet",
     "8 -> 3 3bar",
     OctetTripletTriplet);
  static SwitchOption interfaceColourStructureTripletOctetTriplet
    (interfaceColourStructure,
     "TripletOctetTriplet",
     "3 -> 8 3",
     TripletOctetTriplet);
  static SwitchOption interfaceColourStructureChargedChargedNeutral
    (interfaceColourStructure,
     "ChargedChargedNeutral",
     "q -> q 0",
     ChargedChargedNeutral);
  static SwitchOption interfaceColourStructureNeutralChargedCharged
    (interfaceColourStructure,
     "NeutralChargedCharged",
     "0 -> q qbar",
     NeutralChargedCharged);
  static SwitchOption interfaceColourStructureChargedNeutralCharged
    (interfaceColourStructure,
     "ChargedNeutralCharged",
     "q -> 0 q",
     ChargedNeutralCharged);

  static Switch<SplittingFunction,ShowerInteraction::Type> 
    interfaceInteractionType
    ("InteractionType",
     "Type of the interaction",
     &SplittingFunction::_interactionType, 
     ShowerInteraction::UNDEFINED, false, false);
  static SwitchOption interfaceInteractionTypeQCD
    (interfaceInteractionType,
     "QCD","QCD",ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionTypeQED
    (interfaceInteractionType,
     "QED","QED",ShowerInteraction::QED);

}

void SplittingFunction::persistentOutput(PersistentOStream & os) const {
  using namespace ShowerInteraction;
   os << oenum(_interactionType) << _interactionorder 
      << oenum(_colourStructure) << _colourFactor;
}

void SplittingFunction::persistentInput(PersistentIStream & is, int) {
  using namespace ShowerInteraction;
  is >> ienum(_interactionType) >> _interactionorder 
     >>	ienum(_colourStructure) >> _colourFactor;
}

void SplittingFunction::colourConnection(tShowerParticlePtr parent,
					 tShowerParticlePtr first,
					 tShowerParticlePtr second,
					 const bool back) const {
  if(_colourStructure==TripletTripletOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // ensure input consistency
      assert((!cparent.first &&  cparent.second) || 
	     ( cparent.first && !cparent.second));
      // q -> q g
      if(cparent.first) {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(second);
	newline->addColoured     ( first);
	newline->addAntiColoured (second);
      }
      // qbar -> qbar g
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.second->addAntiColoured(second);
	newline->addColoured(second);
	newline->addAntiColoured(first);
      }
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // ensure input consistency
      assert(( cfirst.first && !cfirst.second) ||
	     (!cfirst.first &&  cfirst.second)); 
      // q -> q g
      if(cfirst.first) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addAntiColoured(second);
	newline->addColoured(second);
	newline->addColoured(parent);
      }
      // qbar -> qbar g
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.second->addColoured(second);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
    }
  }
  else if(_colourStructure==OctetOctetOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      // Randomly decide which of the two gluon products take the
      // colour line passing for the colour of the parent gluon
      // (the other will take the one passing for the anticolour of
      //  the parent gluon).
      if(UseRandom::rndbool()) {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(first);
	cparent.second->addAntiColoured(second);
	newline->addColoured(second);
	newline->addAntiColoured(first);
      }
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(second);
	cparent.second->addAntiColoured(first);
	newline->addColoured(first);
	newline->addAntiColoured(second);
      }
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // Randomly decide which of the two gluon products take the
      // colour line passing for the colour of the parent gluon
      // (the other will take the one passing for the anticolour of
      //  the parent gluon).
      if (UseRandom::rndbool()) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addColoured(parent);
	cfirst.second->addColoured(second);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addAntiColoured(second);
	cfirst.second->addAntiColoured(parent);
	newline->addColoured(parent);
	newline->addColoured(second);
      }
    }
  }
  else if(_colourStructure == OctetTripletTriplet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      cparent.first ->addColoured    ( first);
      cparent.second->addAntiColoured(second);
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // ensure input consistency
      assert(( cfirst.first && !cfirst.second) ||
	     (!cfirst.first &&  cfirst.second));
      // g -> q qbar
      if(cfirst.first) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addColoured(parent);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
      // g -> qbar q
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.second->addAntiColoured(parent);
	newline->addColoured(second);
	newline->addColoured(parent);
      }
    }
  }
  else if(_colourStructure == TripletOctetTriplet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // ensure input consistency
      assert(( cparent.first && !cparent.second) || 
	     (!cparent.first &&  cparent.second));
      // q -> g q
      if(cparent.first) {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(first);
	newline->addColoured    (second);
	newline->addAntiColoured( first);
      }
      // qbar -> g qbar
      else {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.second->addAntiColoured(first);
	newline->addColoured    ( first);
	newline->addAntiColoured(second);
      }
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // q -> g q
      if(parent->id()>0) {
	cfirst.first ->addColoured(parent);
	cfirst.second->addColoured(second);
      }
      else {
	cfirst.first ->addAntiColoured(second);
	cfirst.second->addAntiColoured(parent);
      }
    }
  }
  else if(_colourStructure == ChargedChargedNeutral) {
    if(!parent->data().coloured()) return;
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // q -> q g
      if(cparent.first) {
	cparent.first->addColoured(first);
      }
      // qbar -> qbar g
      if(cparent.second) {
	cparent.second->addAntiColoured(first);
      }
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // q -> q g
      if(cfirst.first) {
	cfirst.first->addColoured(parent);
      }
      // qbar -> qbar g
      if(cfirst.second) {
	cfirst.second->addAntiColoured(parent);
      }
    }
  }
  else if(_colourStructure == ChargedNeutralCharged) {
    if(!parent->data().coloured()) return;
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
				      parent->antiColourLine());
      // q -> q g
      if(cparent.first) {
	cparent.first->addColoured(second);
      }
      // qbar -> qbar g
      if(cparent.second) {
	cparent.second->addAntiColoured(second);
      }
    }
    else {
      ColinePair csecond = ColinePair(second->colourLine(), 
				      second->antiColourLine());
      // q -> q g
      if(csecond.first) {
	csecond.first->addColoured(parent);
      }
      // qbar -> qbar g
      else {
	csecond.second->addAntiColoured(parent);
      }
    }
  }
  else {
    assert(false);
  }
}

void SplittingFunction::doinit() {
  Interfaced::doinit();
  assert(_interactionType!=ShowerInteraction::UNDEFINED);
  assert((_colourStructure>0&&_interactionType==ShowerInteraction::QCD) ||
	 (_colourStructure<0&&_interactionType==ShowerInteraction::QED) );
  if(_colourFactor>0.) return;
  // compute the colour factors if need
  if(_colourStructure==TripletTripletOctet) {
    _colourFactor = 4./3.;
  }
  else if(_colourStructure==OctetOctetOctet) {
    _colourFactor = 3.;
  }
  else if(_colourStructure==OctetTripletTriplet) {
    _colourFactor = 0.5;
  }
  else if(_colourStructure==TripletOctetTriplet) {
    _colourFactor = 4./3.;
  }
  else if(_colourStructure<0) {
    _colourFactor = 1.;
  }
  else {
    assert(false);
  }
}

bool SplittingFunction::checkColours(const IdList & ids) const {
  tcPDPtr pd[3]={getParticleData(ids[0]),
		 getParticleData(ids[1]),
		 getParticleData(ids[2])};
  if(_colourStructure==TripletTripletOctet) {
    if(ids[0]!=ids[1]) return false;
    if((pd[0]->iColour()==PDT::Colour3||pd[0]->iColour()==PDT::Colour3bar) &&
       pd[2]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(_colourStructure==OctetOctetOctet) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(pd[ix]->iColour()!=PDT::Colour8) return false;
    }
    return true;
  }
  else if(_colourStructure==OctetTripletTriplet) {
    if(pd[0]->iColour()!=PDT::Colour8) return false;
    if(pd[1]->iColour()==PDT::Colour3&&pd[2]->iColour()==PDT::Colour3bar)
      return true;
    if(pd[1]->iColour()==PDT::Colour3bar&&pd[2]->iColour()==PDT::Colour3)
      return true;
    return false;
  }
  else if(_colourStructure==TripletOctetTriplet) {
    if(ids[0]!=ids[2]) return false;
    if((pd[0]->iColour()==PDT::Colour3||pd[0]->iColour()==PDT::Colour3bar) &&
       pd[1]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(_colourStructure==ChargedChargedNeutral) {
    if(ids[0]!=ids[1]) return false;
    if(pd[2]->iCharge()!=0) return false;
    if(pd[0]->iCharge()==pd[1]->iCharge()) return true;
    return false;
  }
  else if(_colourStructure==ChargedNeutralCharged) {
    if(ids[0]!=ids[2]) return false;
    if(pd[1]->iCharge()!=0) return false;
    if(pd[0]->iCharge()==pd[2]->iCharge()) return true;
    return false;
  }
  else if(_colourStructure==NeutralChargedCharged) {
    if(ids[1]!=-ids[2]) return false;
    if(pd[0]->iCharge()!=0) return false;
    if(pd[1]->iCharge()==-pd[2]->iCharge()) return true;
    return false;
  }
  else {
    assert(false);
  }
  return false;
}
