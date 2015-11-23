// -*- C++ -*-
//
// MatchboxAmplitudehqqbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehqqbarg class.
//

#include "MatchboxAmplitudehqqbarg.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehqqbarg::MatchboxAmplitudehqqbarg() 
  : interfaceTHooft(126*GeV) {}

MatchboxAmplitudehqqbarg::~MatchboxAmplitudehqqbarg() {}

IBPtr MatchboxAmplitudehqqbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehqqbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehqqbarg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);;
}

void MatchboxAmplitudehqqbarg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehqqbarg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator q=xproc.begin();
  for (; q!=xproc.end(); ++q){
    if ((**q).id()==1 || 
        (**q).id()==2 ||
        (**q).id()==3 ||
        (**q).id()==4 ||
        (**q).id()==5) 
    break;
  }
  if(q==xproc.end()) return false;
  int qid = (**q).id(); 
  xproc.erase(q);
  PDVector::iterator qb=xproc.begin();
  for (; qb!=xproc.end(); ++qb){
    if ((**qb).id()==-qid ){break;}
  }
  if(qb==xproc.end()){ return false;}
  xproc.erase(qb);
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21){xproc.erase(g); break; }
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehqqbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));
 
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehqqbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  unsigned int g=0;
  cPDVector x=amplitudePartonData();
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;}
  for (;g<amplitudePartonData().size();++g){if (x[g]->id() == 21) break;}
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  double alphaS = SM().alphaS();
  double v= 2*MW/gw/sqrt(lastSHat()) ;
  double c = gs*alphaS/3/sqrt(2.)/Constants::pi/v;
  
  if(hel[qbar]==hel[q]){
    largeN = 0;
    return(largeN);
  }
  
  if(hel[qbar]==+1 && hel[q]==-1 && hel[g]==-1){
    largeN = -c*plusProduct(qbar,g)*plusProduct(qbar,g)/(plusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==1 && hel[q]==-1 && hel[g]==1){
    largeN = -c*minusProduct(q,g)*minusProduct(q,g)/(minusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==-1){
    largeN = c*plusProduct(q,g)*plusProduct(q,g)/(plusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==1){
    largeN = c*minusProduct(qbar,g)*minusProduct(qbar,g)/(minusProduct(q,qbar));
    return(largeN);
  }
  // Unknown helicity configuration
  assert(false);
  return(0.);
}

/*Complex MatchboxAmplitudehqqbarg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehqqbarg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehqqbarg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceTHooft,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehqqbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehqqbarg("Herwig::MatchboxAmplitudehqqbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehqqbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudehqqbarg> documentation
    ("MatchboxAmplitudehqqbarg");
  static Parameter<MatchboxAmplitudehqqbarg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehqqbarg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
}

