// -*- C++ -*-
//
// SplittingFunction.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingFunction class.
//

#include "SplittingFunction.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractClass<SplittingFunction,Interfaced>
describeSplittingFunction ("Herwig::SplittingFunction","");

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("The SplittingFunction class is the based class for 1->2 splitting functions"
     " in Herwig");

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
  static SwitchOption interfaceColourStructureSextetSextetOctet
    (interfaceColourStructure,
     "SextetSextetOctet",
     "6 -> 6 8",
     SextetSextetOctet);
  
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

  static Switch<SplittingFunction,ShowerInteraction> 
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

  static Switch<SplittingFunction,bool> interfaceAngularOrdered
    ("AngularOrdered",
     "Whether or not this interaction is angular ordered, "
     "normally only g->q qbar and gamma-> f fbar are the only ones which aren't.",
     &SplittingFunction::angularOrdered_, true, false, false);
  static SwitchOption interfaceAngularOrderedYes
    (interfaceAngularOrdered,
     "Yes",
     "Interaction is angular ordered",
     true);
  static SwitchOption interfaceAngularOrderedNo
    (interfaceAngularOrdered,
     "No",
     "Interaction isn't angular ordered",
     false);

  static Switch<SplittingFunction,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "The scale choice to be used",
     &SplittingFunction::scaleChoice_, 2, false, false);
  static SwitchOption interfaceScaleChoicepT
    (interfaceScaleChoice,
     "pT",
     "pT of the branching",
     0);
  static SwitchOption interfaceScaleChoiceQ2
    (interfaceScaleChoice,
     "Q2",
     "Q2 of the branching",
     1);
  static SwitchOption interfaceScaleChoiceFromAngularOrdering
    (interfaceScaleChoice,
     "FromAngularOrdering",
     "If angular order use pT, otherwise Q2",
     2);

}

void SplittingFunction::persistentOutput(PersistentOStream & os) const {
   os << oenum(_interactionType) << _interactionOrder 
      << oenum(_colourStructure) << _colourFactor
      << angularOrdered_ << scaleChoice_;
}

void SplittingFunction::persistentInput(PersistentIStream & is, int) {
  is >> ienum(_interactionType) >> _interactionOrder 
     >>	ienum(_colourStructure) >> _colourFactor
     >> angularOrdered_ >> scaleChoice_;
}

void SplittingFunction::colourConnection(tShowerParticlePtr parent,
                                         tShowerParticlePtr first,
                                         tShowerParticlePtr second,
					 ShowerPartnerType partnerType, 
                                         const bool back) const {
  if(_colourStructure==TripletTripletOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert((  cparent.first && !cparent.second && 
		partnerType==ShowerPartnerType::QCDColourLine) || 
             ( !cparent.first &&  cparent.second && 
		partnerType==ShowerPartnerType::QCDAntiColourLine));
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
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert((  cfirst.first && !cfirst.second && 
		partnerType==ShowerPartnerType::QCDColourLine) || 
             ( !cfirst.first &&  cfirst.second && 
		partnerType==ShowerPartnerType::QCDAntiColourLine));
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
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
    }
  }
  else if(_colourStructure==OctetOctetOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      // ensure first gluon is hardest
      if( first->id()==second->id() && parent->showerKinematics()->z()<0.5 ) 
	swap(first,second);
      // colour line radiates
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	// The colour line is radiating
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(second);
	cparent.second->addAntiColoured(first);
	newline->addColoured(first);
	newline->addAntiColoured(second);
      }
      // anti colour line radiates
      else if(partnerType==ShowerPartnerType::QCDAntiColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cparent.first->addColoured(first);
	cparent.second->addAntiColoured(second);
	newline->addColoured(second);
	newline->addAntiColoured(first);
      }
      else
	assert(false);
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // The colour line is radiating
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addAntiColoured(second);
	cfirst.second->addAntiColoured(parent);
	newline->addColoured(parent);
	newline->addColoured(second);
      }
      // anti colour line radiates
      else if(partnerType==ShowerPartnerType::QCDAntiColourLine) {
	ColinePtr newline=new_ptr(ColourLine());
	cfirst.first->addColoured(parent);
	cfirst.second->addColoured(second);
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
      else
	assert(false);
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
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
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
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
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
      // Set progenitor
      first->progenitor(parent->progenitor());
      second->progenitor(parent->progenitor());
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
      // Set progenitor
      parent->progenitor(first->progenitor());
      second->progenitor(first->progenitor());
    }
  }
  else if(_colourStructure==SextetSextetOctet) {
    //make sure we're not doing backward evolution
    assert(!back);

    //make sure something sensible
    assert(parent->colourLine() || parent->antiColourLine());
   
    //get the colour lines or anti-colour lines
    bool isAntiColour=true;
    ColinePair cparent;
    if(parent->colourLine()) {
      cparent = ColinePair(const_ptr_cast<tColinePtr>(parent->colourInfo()->colourLines()[0]), 
			   const_ptr_cast<tColinePtr>(parent->colourInfo()->colourLines()[1]));
      isAntiColour=false;
    }
    else {
      cparent = ColinePair(const_ptr_cast<tColinePtr>(parent->colourInfo()->antiColourLines()[0]), 
			   const_ptr_cast<tColinePtr>(parent->colourInfo()->antiColourLines()[1]));
    }
    
    //check for sensible input
    //    assert(cparent.first && cparent.second);

    // sextet has 2 colour lines
    if(!isAntiColour) {
      //pick at random which of the colour topolgies to take
      double topology = UseRandom::rnd();
      if(topology < 0.25) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(second);
        cparent.second->addColoured(first);
        newline->addColoured(first);
        newline->addAntiColoured(second);
      }
      else if(topology >=0.25 && topology < 0.5) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(first);
        cparent.second->addColoured(second);
        newline->addColoured(first);
        newline->addAntiColoured(second); 
      }
      else if(topology >= 0.5 && topology < 0.75) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(second);
        cparent.second->addColoured(first); 
        newline->addColoured(first); 
        newline->addAntiColoured(second); 
      }
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addColoured(first);
        cparent.second->addColoured(second);
        newline->addColoured(first);
        newline->addAntiColoured(second);
      }
    }
    // sextet has 2 anti-colour lines
    else {
      double topology = UseRandom::rnd();
      if(topology < 0.25){
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(second);
        cparent.second->addAntiColoured(first);
        newline->addAntiColoured(first);
        newline->addColoured(second);
      }
      else if(topology >=0.25 && topology < 0.5) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(first);
        cparent.second->addAntiColoured(second);
        newline->addAntiColoured(first);
        newline->addColoured(second); 
      }
      else if(topology >= 0.5 && topology < 0.75) {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(second);
        cparent.second->addAntiColoured(first);
        newline->addAntiColoured(first);
        newline->addColoured(second); 
      }
      else {
        ColinePtr newline=new_ptr(ColourLine());
        cparent.first->addAntiColoured(first);
        cparent.second->addAntiColoured(second);
        newline->addAntiColoured(first);
        newline->addColoured(second);
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
      if (second->dataPtr()->iColour()==PDT::Colour3 ) {
	ColinePtr newline=new_ptr(ColourLine());
	newline->addColoured(second);
	newline->addColoured(parent);
      }
      else if (second->dataPtr()->iColour()==PDT::Colour3bar ) {
	ColinePtr newline=new_ptr(ColourLine());
	newline->addAntiColoured(second);
	newline->addAntiColoured(parent);
      }
    }
  }
  else if(_colourStructure == NeutralChargedCharged ) {
    if(!back) {
      if(first->dataPtr()->coloured()) {
	ColinePtr newline=new_ptr(ColourLine());
	if(first->dataPtr()->iColour()==PDT::Colour3) {
	  newline->addColoured    (first );
	  newline->addAntiColoured(second);
	}
	else if (first->dataPtr()->iColour()==PDT::Colour3bar) {
	  newline->addColoured    (second);
	  newline->addAntiColoured(first );
	}
	else
	  assert(false);
      }
    }
    else {   
      ColinePair cfirst = ColinePair(first->colourLine(), 
				     first->antiColourLine());
      // gamma -> q qbar
      if(cfirst.first) {
	cfirst.first->addAntiColoured(second);
      }
      // gamma -> qbar q
      else if(cfirst.second) {
	cfirst.second->addColoured(second);
      }
      else 
	assert(false);
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
  else if(_colourStructure==SextetSextetOctet) {
    _colourFactor = 10./3.;
  }
  else if(_colourStructure<0) {
    _colourFactor = 1.;
  }
  else {
    assert(false);
  }
}

bool SplittingFunction::checkColours(const IdList & ids) const {
  if(_colourStructure==TripletTripletOctet) {
    if(ids[0]!=ids[1]) return false;
    if((ids[0]->iColour()==PDT::Colour3||ids[0]->iColour()==PDT::Colour3bar) &&
       ids[2]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(_colourStructure==OctetOctetOctet) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ids[ix]->iColour()!=PDT::Colour8) return false;
    }
    return true;
  }
  else if(_colourStructure==OctetTripletTriplet) {
    if(ids[0]->iColour()!=PDT::Colour8) return false;
    if(ids[1]->iColour()==PDT::Colour3&&ids[2]->iColour()==PDT::Colour3bar)
      return true;
    if(ids[1]->iColour()==PDT::Colour3bar&&ids[2]->iColour()==PDT::Colour3)
      return true;
    return false;
  }
  else if(_colourStructure==TripletOctetTriplet) {
    if(ids[0]!=ids[2]) return false;
    if((ids[0]->iColour()==PDT::Colour3||ids[0]->iColour()==PDT::Colour3bar) &&
       ids[1]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(_colourStructure==SextetSextetOctet) {
    if(ids[0]!=ids[1]) return false;
    if((ids[0]->iColour()==PDT::Colour6 || ids[0]->iColour()==PDT::Colour6bar) &&
       ids[2]->iColour()==PDT::Colour8) return true;
    return false;
  }
  else if(_colourStructure==ChargedChargedNeutral) {
    if(ids[0]!=ids[1]) return false;
    if(ids[2]->iCharge()!=0) return false;
    if(ids[0]->iCharge()==ids[1]->iCharge()) return true;
    return false;
  }
  else if(_colourStructure==ChargedNeutralCharged) {
    if(ids[0]!=ids[2]) return false;
    if(ids[1]->iCharge()!=0) return false;
    if(ids[0]->iCharge()==ids[2]->iCharge()) return true;
    return false;
  }
  else if(_colourStructure==NeutralChargedCharged) {
    if(ids[1]->id()!=-ids[2]->id()) return false;
    if(ids[0]->iCharge()!=0) return false;
    if(ids[1]->iCharge()==-ids[2]->iCharge()) return true;
    return false;
  }
  else {
    assert(false);
  }
  return false;
}

namespace {

  bool hasColour(tPPtr p) {
    PDT::Colour colour = p->dataPtr()->iColour();
    return colour==PDT::Colour3 || colour==PDT::Colour8 || colour == PDT::Colour6;
  }

  bool hasAntiColour(tPPtr p) {
    PDT::Colour colour = p->dataPtr()->iColour();
    return colour==PDT::Colour3bar || colour==PDT::Colour8 || colour == PDT::Colour6bar;
  }
  
}

void SplittingFunction::evaluateFinalStateScales(ShowerPartnerType partnerType,
						 Energy scale, double z,
						 tShowerParticlePtr parent,
						 tShowerParticlePtr emitter,
						 tShowerParticlePtr emitted) {
  // identify emitter and emitted
  double zEmitter = z, zEmitted = 1.-z;
  bool bosonSplitting(false);
  // special for g -> gg, particle highest z is emitter
  if(emitter->id() == emitted->id() && emitter->id() == parent->id() &&
     zEmitted > zEmitter) {
    swap(zEmitted,zEmitter);
    swap( emitted, emitter);
  }
  // otherwise if particle ID same
  else if(emitted->id()==parent->id()) {
    swap(zEmitted,zEmitter);
    swap( emitted, emitter);
  }
  // no real emitter/emitted
  else if(emitter->id()!=parent->id()) {
    bosonSplitting = true;
  }
  // may need to add angularOrder flag here
  // now the various scales
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    assert(colourStructure()==ChargedChargedNeutral ||
	   colourStructure()==ChargedNeutralCharged ||
	   colourStructure()==NeutralChargedCharged );
    // normal case
    if(!bosonSplitting) {
      assert(colourStructure()==ChargedChargedNeutral ||
	     colourStructure()==ChargedNeutralCharged );
      // set the scales
      // emitter
      emitter->scales().QED         = zEmitter*scale;
      emitter->scales().QED_noAO    =          scale;
      emitter->scales().QCD_c       = min(scale,parent->scales().QCD_c      );
      emitter->scales().QCD_c_noAO  = min(scale,parent->scales().QCD_c_noAO );
      emitter->scales().QCD_ac      = min(scale,parent->scales().QCD_ac     );
      emitter->scales().QCD_ac_noAO = min(scale,parent->scales().QCD_ac_noAO);
      // emitted 
      emitted->scales().QED         = zEmitted*scale;
      emitted->scales().QED_noAO    =          scale;
      emitted->scales().QCD_c       = ZERO;
      emitted->scales().QCD_c_noAO  = ZERO;
      emitted->scales().QCD_ac      = ZERO;
      emitted->scales().QCD_ac_noAO = ZERO;
    }
    // gamma -> f fbar
    else {
      assert(colourStructure()==NeutralChargedCharged );
      // emitter
      emitter->scales().QED         = zEmitter*scale;
      emitter->scales().QED_noAO    =          scale;
      if(hasColour(emitter)) {
	emitter->scales().QCD_c       = zEmitter*scale;
	emitter->scales().QCD_c_noAO  =          scale;
      }
      if(hasAntiColour(emitter)) {
	emitter->scales().QCD_ac      = zEmitter*scale;
	emitter->scales().QCD_ac_noAO =          scale;
      }
      // emitted 
      emitted->scales().QED         = zEmitted*scale;
      emitted->scales().QED_noAO    =          scale;
      if(hasColour(emitted)) {
	emitted->scales().QCD_c       = zEmitted*scale;
	emitted->scales().QCD_c_noAO  =          scale;
      }
      if(hasAntiColour(emitted)) {
	emitted->scales().QCD_ac      = zEmitted*scale;
	emitted->scales().QCD_ac_noAO =          scale;
      }
    }
  }
  // QCD
  else {
   // normal case eg q -> q g and g -> g g
    if(!bosonSplitting) {
      emitter->scales().QED         = min(scale,parent->scales().QED     );
      emitter->scales().QED_noAO    = min(scale,parent->scales().QED_noAO);
      if(partnerType==ShowerPartnerType::QCDColourLine) {
	emitter->scales().QCD_c       = zEmitter*scale;
	emitter->scales().QCD_c_noAO  =          scale;
	emitter->scales().QCD_ac      = min(zEmitter*scale,parent->scales().QCD_ac     );
	emitter->scales().QCD_ac_noAO = min(         scale,parent->scales().QCD_ac_noAO);
      }
      else {
	emitter->scales().QCD_c       = min(zEmitter*scale,parent->scales().QCD_c      );
	emitter->scales().QCD_c_noAO  = min(         scale,parent->scales().QCD_c_noAO );
	emitter->scales().QCD_ac      = zEmitter*scale;
	emitter->scales().QCD_ac_noAO =          scale;
      }
      // emitted 
      emitted->scales().QED         = ZERO;
      emitted->scales().QED_noAO    = ZERO;
      emitted->scales().QCD_c       = zEmitted*scale;
      emitted->scales().QCD_c_noAO  =          scale;
      emitted->scales().QCD_ac      = zEmitted*scale;
      emitted->scales().QCD_ac_noAO =          scale;
    }
    // g -> q qbar
    else {
      // emitter
      if(emitter->dataPtr()->charged()) {
	emitter->scales().QED         = zEmitter*scale;
	emitter->scales().QED_noAO    =          scale;
      }
      emitter->scales().QCD_c       = zEmitter*scale;
      emitter->scales().QCD_c_noAO  =          scale;
      emitter->scales().QCD_ac      = zEmitter*scale;
      emitter->scales().QCD_ac_noAO =          scale;
      // emitted 
      if(emitted->dataPtr()->charged()) {
	emitted->scales().QED         = zEmitted*scale;
	emitted->scales().QED_noAO    =          scale;
      }
      emitted->scales().QCD_c       = zEmitted*scale;
      emitted->scales().QCD_c_noAO  =          scale;
      emitted->scales().QCD_ac      = zEmitted*scale;
      emitted->scales().QCD_ac_noAO =          scale;
    }
  }
}

void SplittingFunction::evaluateInitialStateScales(ShowerPartnerType partnerType,
						   Energy scale, double z,
						   tShowerParticlePtr parent,
						   tShowerParticlePtr spacelike,
						   tShowerParticlePtr timelike) {
  // scale for time-like child
  Energy AOScale = (1.-z)*scale;
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    if(parent->id()==spacelike->id()) {
      // parent
      parent   ->scales().QED         =   scale;
      parent   ->scales().QED_noAO    =   scale;
      parent   ->scales().QCD_c       = min(scale,spacelike->scales().QCD_c      );
      parent   ->scales().QCD_c_noAO  = min(scale,spacelike->scales().QCD_c_noAO );
      parent   ->scales().QCD_ac      = min(scale,spacelike->scales().QCD_ac     );
      parent   ->scales().QCD_ac_noAO = min(scale,spacelike->scales().QCD_ac_noAO);
      // timelike
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      timelike->scales().QCD_c       =    ZERO;
      timelike->scales().QCD_c_noAO  =    ZERO;
      timelike->scales().QCD_ac      =    ZERO;
      timelike->scales().QCD_ac_noAO =    ZERO;
    }
    else if(parent->id()==timelike->id()) {
      parent   ->scales().QED         =   scale;
      parent   ->scales().QED_noAO    =   scale;
      if(hasColour(parent)) {
	parent   ->scales().QCD_c       = scale;
	parent   ->scales().QCD_c_noAO  = scale;
      }
      if(hasAntiColour(parent)) {
	parent   ->scales().QCD_ac      = scale;
	parent   ->scales().QCD_ac_noAO = scale;
      }
      // timelike 
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      if(hasColour(timelike)) {
	timelike->scales().QCD_c       = AOScale;
	timelike->scales().QCD_c_noAO  =   scale;
      }
      if(hasAntiColour(timelike)) {
	timelike->scales().QCD_ac      = AOScale;
	timelike->scales().QCD_ac_noAO =   scale;
      }
    }
    else {
      parent   ->scales().QED         = scale;
      parent   ->scales().QED_noAO    = scale;
      parent   ->scales().QCD_c       = ZERO ;
      parent   ->scales().QCD_c_noAO  = ZERO ;
      parent   ->scales().QCD_ac      = ZERO ;
      parent   ->scales().QCD_ac_noAO = ZERO ;
      // timelike 
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
      if(hasColour(timelike)) {
	timelike->scales().QCD_c       = min(AOScale,spacelike->scales().QCD_ac     );
	timelike->scales().QCD_c_noAO  = min(  scale,spacelike->scales().QCD_ac_noAO);
      }
      if(hasAntiColour(timelike)) {
	timelike->scales().QCD_ac      = min(AOScale,spacelike->scales().QCD_c      );
	timelike->scales().QCD_ac_noAO = min(  scale,spacelike->scales().QCD_c_noAO );
      }
    }
  }
  // QCD
  else {
    // timelike 
    if(timelike->dataPtr()->charged()) {
      timelike->scales().QED         = AOScale;
      timelike->scales().QED_noAO    =   scale;
    }
    if(hasColour(timelike)) {
      timelike->scales().QCD_c       = AOScale;
      timelike->scales().QCD_c_noAO  =   scale;
    }
    if(hasAntiColour(timelike)) {
      timelike->scales().QCD_ac      = AOScale;
      timelike->scales().QCD_ac_noAO =   scale;
    }
    if(parent->id()==spacelike->id()) {
      parent   ->scales().QED         = min(scale,spacelike->scales().QED        );
      parent   ->scales().QED_noAO    = min(scale,spacelike->scales().QED_noAO   );
      parent   ->scales().QCD_c       = min(scale,spacelike->scales().QCD_c      );
      parent   ->scales().QCD_c_noAO  = min(scale,spacelike->scales().QCD_c_noAO );
      parent   ->scales().QCD_ac      = min(scale,spacelike->scales().QCD_ac     );
      parent   ->scales().QCD_ac_noAO = min(scale,spacelike->scales().QCD_ac_noAO);
    }
    else {
      if(parent->dataPtr()->charged()) {
	parent   ->scales().QED         = scale;
	parent   ->scales().QED_noAO    = scale;
      }
      if(hasColour(parent)) {
	parent   ->scales().QCD_c      = scale;
	parent   ->scales().QCD_c_noAO  = scale;
      }
      if(hasAntiColour(parent)) {
	parent   ->scales().QCD_ac      = scale;
	parent   ->scales().QCD_ac_noAO = scale;
      }
    }
  }
}

void SplittingFunction::evaluateDecayScales(ShowerPartnerType partnerType,
					    Energy scale, double z,
					    tShowerParticlePtr parent,
					    tShowerParticlePtr spacelike,
					    tShowerParticlePtr timelike) {
  assert(parent->id()==spacelike->id());
  // angular-ordered scale for 2nd child
  Energy AOScale = (1.-z)*scale;
  // QED
  if(partnerType==ShowerPartnerType::QED) {
    // timelike
    timelike->scales().QED         = AOScale;
    timelike->scales().QED_noAO    =   scale;
    timelike->scales().QCD_c       =    ZERO;
    timelike->scales().QCD_c_noAO  =    ZERO;
    timelike->scales().QCD_ac      =    ZERO;
    timelike->scales().QCD_ac_noAO =    ZERO;
    // spacelike
    spacelike->scales().QED         =   scale;
    spacelike->scales().QED_noAO    =   scale;
  }
  // QCD
  else {
    // timelike 
    timelike->scales().QED         = ZERO;
    timelike->scales().QED_noAO    = ZERO;
    timelike->scales().QCD_c       = AOScale;
    timelike->scales().QCD_c_noAO  =   scale;
    timelike->scales().QCD_ac      = AOScale;
    timelike->scales().QCD_ac_noAO =   scale;
    // spacelike
    spacelike->scales().QED         = max(scale,parent->scales().QED        );
    spacelike->scales().QED_noAO    = max(scale,parent->scales().QED_noAO   );
  }
  spacelike->scales().QCD_c       = max(scale,parent->scales().QCD_c      );
  spacelike->scales().QCD_c_noAO  = max(scale,parent->scales().QCD_c_noAO );
  spacelike->scales().QCD_ac      = max(scale,parent->scales().QCD_ac     );
  spacelike->scales().QCD_ac_noAO = max(scale,parent->scales().QCD_ac_noAO);
}
