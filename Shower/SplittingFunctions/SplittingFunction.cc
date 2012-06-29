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
  static SwitchOption interfaceColourStructureSextetSextetOctet
    (interfaceColourStructure,
     "SextetSextetOctet",
     "6 -> 6 8",
     SextetSextetOctet);
  

  static Switch<SplittingFunction,ShowerInteraction::Type> 
    interfaceInteractionType
    ("InteractionType",
     "Type of the interaction",
     &SplittingFunction::_interactionType, 
     ShowerInteraction::UNDEFINED, false, false);
  static SwitchOption interfaceInteractionTypeQCD
    (interfaceInteractionType,
     "QCD","QCD",ShowerInteraction::QCD);


  static Switch<SplittingFunction,int> interfaceSplittingColourMethod
    ("SplittingColourMethod",
     "Choice of assigning colour in 8->88 splittings.",
     &SplittingFunction::_splittingColourMethod, 0, false, false);
  static SwitchOption interfaceSplittingColourMethodRandom
    (interfaceSplittingColourMethod,
     "Random",
     "Choose colour assignments randomly.",
     0);
  static SwitchOption interfaceSplittingColourMethodCorrectLines
    (interfaceSplittingColourMethod,
     "CorrectLines",
     "Choose correct lines for colour.",
     1);
  static SwitchOption interfaceSplittingColourMethodRandomRecord
    (interfaceSplittingColourMethod,
     "RandomRecord",
     "Choose colour assignments randomly and record the result.",
     2);
}

void SplittingFunction::persistentOutput(PersistentOStream & os) const {
  using namespace ShowerInteraction;
   os << oenum(_interactionType) << _interactionorder 
      << oenum(_colourStructure) << _colourFactor
      << _splittingColourMethod;
}

void SplittingFunction::persistentInput(PersistentIStream & is, int) {
  using namespace ShowerInteraction;
  is >> ienum(_interactionType) >> _interactionorder 
     >>	ienum(_colourStructure) >> _colourFactor
     >> _splittingColourMethod;
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
      // Set progenitor
      first->setProgenitor(parent->progenitor());
      second->setProgenitor(parent->progenitor());
      // Random radiation choice
      first->setRadiationLine(0);
      second->setRadiationLine(0);      
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
      // Set progenitor
      parent->setProgenitor(first->progenitor());
      second->setProgenitor(first->progenitor());
      // Random radiation choice
      parent->setRadiationLine(0);
      second->setRadiationLine(0); 
    }
  }
  else if(_colourStructure==OctetOctetOctet) {
    if(!back) {
      ColinePair cparent = ColinePair(parent->colourLine(), 
                                      parent->antiColourLine());
      // ensure input consistency
      assert(cparent.first&&cparent.second);
      // The choice of colour line is determined by the
      // radiation line of the parent. 
      // If the radiation line is non-zero and the
      // scale of the parent is above the second scale of the 
      // progenitor it will only radiate from the chosen radiation
      // line. Otherwise the parent will radiate randomly.
      // Initializing radiation lines
      first->setRadiationLine(0);
      second->setRadiationLine(0);
      // Switch to choose random or non-random choice of lines
      bool randomchoice = 0;
      // Radiation line
      int radiationLine = 0;
      if (_splittingColourMethod == 1){
        // Choose the appropriate colour lines
        if ((parent->radiationLine() == 1 || parent->radiationLine() == 2)  && parent->progenitor() ) {
          if (parent->evolutionScale() > parent->progenitor()->evolutionScale2()){
	    // Parent has a radiation line, so the line which should
	    // radiate, and therefore the choice of which colour line
	    // to pass onto which child, is already determined.
	    randomchoice = 1;
            if(parent->radiationLine() == 2) {
              // The anti-colour line is radiating
      	      ColinePtr newline=new_ptr(ColourLine());
      	      cparent.first->addColoured(first);
      	      cparent.second->addAntiColoured(second);
      	      newline->addColoured(second);
      	      newline->addAntiColoured(first);
	      // Set the radiation line for the children
	      radiationLine = parent->radiationLine();
            }	
            else {
              // The colour line is radiating
      	      ColinePtr newline=new_ptr(ColourLine());
      	      cparent.first->addColoured(second);
      	      cparent.second->addAntiColoured(first);
      	      newline->addColoured(first);
      	      newline->addAntiColoured(second);
      	      // Set the radiation line for the children
	      radiationLine = parent->radiationLine();
            } 
	  } 	
        }
      }
      if (randomchoice == 0) {
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
	  if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
	    if (parent->radiationLine() == 1 || parent->radiationLine() == 2){
	      // Record which line radiates
	      radiationLine = 2;
	    }
	  }
        }	
        else {
      	  ColinePtr newline=new_ptr(ColourLine());
      	  cparent.first->addColoured(second);
      	  cparent.second->addAntiColoured(first);
      	  newline->addColoured(first);
          newline->addAntiColoured(second);
	  if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
	    if (parent->radiationLine() == 1 || parent->radiationLine() == 2){
	      // Record which line radiates
	      radiationLine = 1;
	    }
	  }
        }    
      }
      if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
	if (parent->radiationLine() == 1 || parent->radiationLine() == 2){
      	  // Set the radiation line for the children
      	  first->setRadiationLine(radiationLine);
	  second->setRadiationLine(0);
	  // Set the progenitors for the children
	  first->setProgenitor(parent->progenitor());
	  second->setProgenitor(parent->progenitor()); 
	}    
      }            
    }
    else {
      ColinePair cfirst = ColinePair(first->colourLine(), 
                                     first->antiColourLine());
      // ensure input consistency
      assert(cfirst.first&&cfirst.second);
      // The choice of colour line is determined by the
      // radiation line of the parent. 
      // If the radiation line is non-zero and the
      // scale of the parent is above the second scale of the 
      // progenitor it will only radiate from the chosen radiation
      // line. Otherwise the parent will radiate randomly.           
      // Initializing radiation lines
      parent->setRadiationLine(0);
      second->setRadiationLine(0);
      // Switch to choose random or non-random choice of lines
      bool randomchoice = 0;
      // Radiation line
      int radiationLine = 0;
      if (_splittingColourMethod == 1){
        // Choose the appropriate colour lines
        if ((first->radiationLine() == 1 || first->radiationLine() == 2) && first->progenitor()) {
          if (first->evolutionScale() > first->progenitor()->evolutionScale2()){
            // Parent has a radiation line, so the line which should
	    // radiate, and therefore the choice of which colour line
	    // to pass onto which child, is already determined.  
	    randomchoice = 1;
    	    if (first->radiationLine() == 2) {
	      // The anti-colour line is radiating
      	      ColinePtr newline=new_ptr(ColourLine());
      	      cfirst.first->addColoured(parent);
      	      cfirst.second->addColoured(second);
      	      newline->addAntiColoured(second);
      	      newline->addAntiColoured(parent);
              // Set the radiation line for the children
      	      radiationLine = first->radiationLine();

            }    
            else {
	      // The colour line is radiating
      	      ColinePtr newline=new_ptr(ColourLine());
      	      cfirst.first->addAntiColoured(second);
      	      cfirst.second->addAntiColoured(parent);
      	      newline->addColoured(parent);
      	      newline->addColoured(second);
	      // Set the radiation line for the children
	      radiationLine = first->radiationLine();
            }    
	  }
        }
      }
      if (randomchoice == 0) {
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
	  if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
	    if (first->radiationLine() == 1 || first->radiationLine() == 2){
	      // Record which line radiates
	      radiationLine = 2;
	    }
	  }
        }   
        else {
       	  ColinePtr newline=new_ptr(ColourLine());
      	  cfirst.first->addAntiColoured(second);
      	  cfirst.second->addAntiColoured(parent);
      	  newline->addColoured(parent);
      	  newline->addColoured(second);
	  if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
	    if (first->radiationLine() == 1 || first->radiationLine() == 2){
	      // Record which line radiates
	      radiationLine = 1;
	    }
	  }
        }         
      }
      if (_splittingColourMethod == 1 || _splittingColourMethod == 2){
        if (first->radiationLine() == 1 || first->radiationLine() == 2){
	  // Set the radiation line for the children
     	  parent->setRadiationLine(radiationLine);
	  second->setRadiationLine(0);
	  // Set the progenitors for the children
	  parent->setProgenitor(first->progenitor());
	  second->setProgenitor(first->progenitor()); 	  
	}   	      
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
      // Set progenitor
      first->setProgenitor(parent->progenitor());
      second->setProgenitor(parent->progenitor());
      // Random radiation choice
      first->setRadiationLine(0);
      second->setRadiationLine(0); 
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
      parent->setProgenitor(first->progenitor());
      second->setProgenitor(first->progenitor());
      // Random radiation choice
      parent->setRadiationLine(0);
      second->setRadiationLine(0); 
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
      first->setProgenitor(parent->progenitor());
      second->setProgenitor(parent->progenitor());
      // Random radiation choice
      first->setRadiationLine(0);
      second->setRadiationLine(0); 
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
      parent->setProgenitor(first->progenitor());
      second->setProgenitor(first->progenitor());
      // Random radiation choice
      parent->setRadiationLine(0);
      second->setRadiationLine(0); 
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
  else {
    assert(false);
  }
}

void SplittingFunction::doinit() {
  Interfaced::doinit();
  assert(_interactionType!=ShowerInteraction::UNDEFINED);
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
  else if(_colourStructure==SextetSextetOctet) {
    if(ids[0]!=ids[1]) return false;
    if((pd[0]->iColour()==PDT::Colour6 || pd[0]->iColour()==PDT::Colour6bar) &&
       pd[2]->iColour()==PDT::Colour8) return true;
    return false;
  }    
  else {
    assert(false);
  }
  return false;
}
