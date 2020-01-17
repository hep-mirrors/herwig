// -*- C++ -*-
//
// MultiWeight.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiWeight class.
//

#include <cmath>
#include <fstream>

#include "MultiWeight.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <boost/algorithm/string.hpp>

using namespace Herwig;

MultiWeight::MultiWeight() {}


IBPtr MultiWeight::clone() const {
  return new_ptr(*this);
}

IBPtr MultiWeight::fullclone() const {
  return new_ptr(*this);
}


//Define the necessary fortran analyse(s) here. Make sure you include the underscore. 


void MultiWeight::doinitrun() {
  AnalysisHandler::doinitrun();    
}

//The following function analyses an event. You may also add a C++ analysis here. 
void MultiWeight::analyze(tEventPtr event, long, int, int) {

  /* variables to be
   * filled
   */
  double real_weight;
  int iapp = 0;

  /* 
   * loop over the optional weights
   * use the "special" weight values 
   * to identify them
   */ 
  for (map<string,double>::const_iterator it= event->optionalWeights().begin(); it!=event->optionalWeights().end(); ++it){
    
    string one = it->first; // contains the weight information
    double two = it->second; // contains the weight

    // print for debugging
    // cout << "it = " << it->first << "->" << it->second << endl;
    
    /* 
     * aMCFast weights
     */ 
    if(it->second == -111) {
      one.erase(one.begin(),one.begin()+7);
      one.erase(one.end()-1,one.end());
      iapp++;
      cout << "aMCFast: " << endl;
      cout << "\tit = " << it->first << endl;
    }
    
    /* 
     * mg5 clustering info
     */ 
    if(it->second == -222) {
      cout << "mg5 clustering: " << endl;
      cout << "\tit = " << it->first << endl;
    }
    /*
     * mg5 scale clustering info
     */
    if(it->second == -333) {
      cout << "mg5 clustering scale: " << endl;
      cout << "\tit = " << it->first << endl;
    }

    /*
     * mg5 scale clustering info
     */
    if(it->second == -999) {
      cout << "FxFx info: " << endl;
      cout << "\tit = " << it->first << endl;
    }
    
    /* if none of the special tags 
     * then it's an optional weight
     */
    if(it->second != -111 && it->second != -222 && it->second != -333 && it->second != -999) {
      cout << "optional weight:" << endl;
      cout << "\tit = " << it->first << "->" << it->second << endl;
    }
    
    /* 
     * the "central" weight for the event
     */
    string central = "central";
    if(one == central) { real_weight = it->second; }
  }

}


NoPIOClassDescription<MultiWeight> MultiWeight::initMultiWeight;
// Definition of the static class description member.

void MultiWeight::Init() {
  std::cout << "MultiWeight analysis loaded" << endl;
 
  static ClassDocumentation<MultiWeight> documentation
    ("The MultiWeight class prints is an interface to a Multiple Weight analysis");

}

void MultiWeight::dofinish() {
  AnalysisHandler::dofinish();
}

