// -*- C++ -*-
//
// MatchboxAmplitudehgg.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehgg class.
//



#include "MatchboxAmplitudehgg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"


using namespace Herwig;

MatchboxAmplitudehgg::MatchboxAmplitudehgg()
  : interfaceTHooft(126*GeV) {}

MatchboxAmplitudehgg::~MatchboxAmplitudehgg() {}

IBPtr MatchboxAmplitudehgg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehgg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehgg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->mass();
  nPoints(3);
}

void MatchboxAmplitudehgg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->mass();
  nPoints(3);
}

bool MatchboxAmplitudehgg::canHandle(const PDVector& proc) const {
  //cout<<"Amplitudehgg::canHandle";
  if ( proc.size() != 3 ) return false;
  PDVector xproc = proc;
  //for (int i=0;i<3;i++){cout<<proc[i]->id();}
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21) {xproc.erase(g); --g;}
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehgg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  //cout<<"prepare erreicht"<<endl;
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));
  //cout<<"momenta werden gesetzt"<<endl;
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehgg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
 //cout<<"Amplitudehgg::evaluate";
 unsigned int g1=0;
 unsigned int g2=0;
 cPDVector x=amplitudePartonData();
 for (;g1<amplitudePartonData().size();++g1){
   if (x[g1]->id()==21) {
     for (g2=g1+1;g2<amplitudePartonData().size();++g2){if (x[g2]->id()==21) break;}
     break;
   }
 }
 //cout<<"g1: "<<g1<<endl;
 //cout<<"g2: "<<g2<<endl;
 if (g1==g2) cout<<"Fehler bei der Teilchenzuordnung!!!"<<endl;
 // double gw = sqrt(4*Constants::pi*0.00775855) / sqrt(SM().sin2ThetaW());
 double gw = sqrt(4*Constants::pi*SM().alphaEM()) / sqrt(SM().sin2ThetaW());
 double v= 2*MW/gw/sqrt(lastSHat()); // Eigentlich invEnergy aber alles auf Shat normiert. -> die amplitude ist double statt Energy!

// CrossingMap Ausgabe:
// for(int i=0; i<3; ++i){cout<<"i="<<i<<":  PartonData[i]: "<<x[i]->id()<<"  Map[i]: "<<lastCrossingMap()[i]<<endl;}       
// cout<<"i=4:  "<<"  Map[i]: "<<lastCrossingMap()[4]<<endl;                                                                
// cout<<endl; 

 double c = SM().alphaS()/3/Constants::pi/v; // nur /3 statt /6 weil vermutlich 1/2 in die farbbasis gezogen wird
// cout<<"hel[g1]"<<hel[g1]<<"  hel[g2]"<<hel[g2]<<endl<<flush;
 if (hel[g1]==-hel[g2]){
     largeN=0;
     return largeN;
 }
 if (hel[g1]==hel[g2] && hel[g1]==1){
     largeN=c*plusProduct(g1,g2)*plusProduct(g1,g2);
//     cout<<"hel[g1]"<<hel[g1]<<"  hel[g2]"<<hel[g2]<<" hgg largeN = "<<largeN<<endl;
     return largeN;
 }
 if (hel[g1]==hel[g2] && hel[g1]==-1){
     largeN=c*minusProduct(g1,g2)*minusProduct(g1,g2);
//     cout<<"hel[g1]"<<hel[g1]<<"  hel[g2]"<<hel[g2]<<" hgg largeN = "<<largeN<<endl;
     return largeN;
 }
 cout<<"unbedachte Helizitaetskonfiguration"<<flush;
 return 0;
}

Complex MatchboxAmplitudehgg::evaluateOneLoop(size_t a , const vector<int>& Hel) {
  Complex E = SM().alphaS()/Constants::pi*11/4;
  Complex largeN = 0.;
  return (E*evaluate(a,Hel,largeN));   
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehgg::persistentOutput(PersistentOStream &os) const {
    os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehgg::persistentInput(PersistentIStream &is, int) {
    is >> iunit(interfaceTHooft,GeV);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehgg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehgg("Herwig::MatchboxAmplitudehgg", "HwMatchbox.so");

void MatchboxAmplitudehgg::Init() {

  static ClassDocumentation<MatchboxAmplitudehgg> documentation
    ("MatchboxAmplitudehgg");
  static Parameter<MatchboxAmplitudehgg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehgg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

