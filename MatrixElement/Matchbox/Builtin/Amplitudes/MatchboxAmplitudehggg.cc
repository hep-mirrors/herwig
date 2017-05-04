// -*- C++ -*-
//
// MatchboxAmplitudehggg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehggg class.
//

#include "MatchboxAmplitudehggg.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehggg::MatchboxAmplitudehggg()
  : interfaceTHooft(126*GeV) {}

MatchboxAmplitudehggg::~MatchboxAmplitudehggg() {}

IBPtr MatchboxAmplitudehggg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehggg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehggg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

void MatchboxAmplitudehggg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehggg::canHandle(const PDVector& proc) const {

  if ( proc.size() != 4 ) return false;
  PDVector xproc = proc;

  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21) {xproc.erase(g); --g;}
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehggg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

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

Complex MatchboxAmplitudehggg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  
  unsigned int p=0;
  unsigned int q=0;
  unsigned int r=0;
  cPDVector x=amplitudePartonData();
  for (;p<amplitudePartonData().size();++p){
    if (x[p]->id()==21) {
      for (q=p+1;q<amplitudePartonData().size();++q){
        if (x[q]->id()==21){
         for (r=q+1;r<amplitudePartonData().size();++r){if (x[r]->id()==21) break;}
          } 
        break;}
      break;
    }
  }
  // Wrong particle assignment. There have to be three distinct Gluons p, q and r.
  assert(!((p==q) || (p==r) || (q==r)));
 
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  double v= 2*MW/gw/sqrt(lastSHat()) ;
  Complex c = Complex (0.,0.);                                                                                                                    
  
  // Assertion makes sure that there is no crossing, which causes a relative minus sign. If Assert -> new, more general code is needed                                                                                                                                    
  assert(amplitudeToColourMap()[1]==0 && amplitudeToColourMap()[2]==1 && amplitudeToColourMap()[3]==2 );

  if (a==0){
    c = Complex (0.,-1.)*sqrt(2.)*SM().alphaS()/3./Constants::pi/v*gs ; 
  } 
  else {
    if (a==1){
      c = Complex (0.,+1.)*sqrt(2.)*SM().alphaS()/3./Constants::pi/v*gs ; 
    }
    else{ 
      //The Colourbasis a is not appropriate for this process. 
      //  hggg ~ f^{abc} -> ~ tr(t^a,t^b,t^c) - tr(t^a,t^c,t^b) -> a in {0,1}
      assert(true);
    }
  }
  if(hel[p]==+1 && hel[q]==+1 && hel[r]==+1){
    largeN = c*(invariant(p,q)*invariant(p,q)+2*invariant(p,q)*invariant(p,r)
            +invariant(p,r)*invariant(p,r)+2*invariant(p,q)*invariant(q,r)
            +2*invariant(p,r)*invariant(q,r)+invariant(q,r)*invariant(q,r))
            /plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==+1 && hel[r]==-1){
    largeN = -c*minusProduct(p,q)*minusProduct(p,q)*minusProduct(p,q)/minusProduct(p,r)/minusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==-1 && hel[r]==+1){
    largeN = -c*minusProduct(p,r)*minusProduct(p,r)*minusProduct(p,r)/minusProduct(p,q)/minusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==-1 && hel[r]==-1){
    largeN = c*plusProduct(q,r)*plusProduct(q,r)*plusProduct(q,r)/plusProduct(p,q)/plusProduct(p,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==+1 && hel[r]==+1){
    largeN = -c*minusProduct(q,r)*minusProduct(q,r)*minusProduct(q,r)/minusProduct(p,q)/minusProduct(p,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==+1 && hel[r]==-1){
    largeN = c*plusProduct(p,r)*plusProduct(p,r)*plusProduct(p,r)/plusProduct(p,q)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==-1 && hel[r]==+1){
    largeN = c*plusProduct(p,q)*plusProduct(p,q)*plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==-1 && hel[r]==-1){
    largeN = -c*(invariant(p,q)*invariant(p,q)+2*invariant(p,q)*invariant(p,r)
            +invariant(p,r)*invariant(p,r)+2*invariant(p,q)*invariant(q,r)
            +2*invariant(p,r)*invariant(q,r)+invariant(q,r)*invariant(q,r))
            /plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  // Unknown helicity configuration
  assert(false);
  return(0.);
}

/*Complex MatchboxAmplitudehggg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehggg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehggg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceTHooft,GeV);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehggg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehggg("Herwig::MatchboxAmplitudehggg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehggg::Init() {

  static ClassDocumentation<MatchboxAmplitudehggg> documentation
    ("MatchboxAmplitudehggg");

  /*  // not used guess leftover from validation (mu2() variation)
  static Parameter<MatchboxAmplitudehggg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehggg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  */
}

 
