// -*- C++ -*-
//
// MatchboxAmplitudehbbbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehbbbar class.
//

#include "MatchboxAmplitudehbbbar.h"
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

MatchboxAmplitudehbbbar::MatchboxAmplitudehbbbar() {}

MatchboxAmplitudehbbbar::~MatchboxAmplitudehbbbar() {}

IBPtr MatchboxAmplitudehbbbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehbbbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehbbbar::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(3);;
}

void MatchboxAmplitudehbbbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(3);
}

bool MatchboxAmplitudehbbbar::canHandle(const PDVector& proc) const {
  if ( proc.size() != 3 )
    return false;
  PDVector xproc = proc;
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
      if ((**qb).id()==-qid) break;
  }
  if(qb==xproc.end()) return false;
  xproc.erase(qb);
  PDVector::iterator h0=xproc.begin();
  if (xproc.size()==1 && (**h0).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehbbbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
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

Complex MatchboxAmplitudehbbbar::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  cPDVector x=amplitudePartonData();   
  /*for (int i=0;i<3;++i){
    cout<<"x["<<i<<"]: "<<x[i]->id()<<endl;
  }*/
  //cout<<"teilchenidentifizierung";
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} //cout<<"x[q]"<<x[q]->id()<<" hel: "<<hel[q]<<endl;
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;} //cout<<"x[qbar]"<<x[qbar]->id()<<" hel: "<<hel[qbar]<<endl;
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  
  long id=x[q]->id();
  Energy Mf = 0*GeV;
  switch(id){
    case 1 : Mf = interfaceDMass; /*cout<<"d"<<ounit(Mf,GeV)<<endl;*/ break;
    case 2 : Mf = interfaceUMass; /*cout<<"u"<<ounit(Mf,GeV)<<endl;*/ break;
    case 3 : Mf = interfaceSMass; /*cout<<"s"<<ounit(Mf,GeV)<<endl;*/ break;
    case 4 : Mf = interfaceCMass; /*cout<<"c"<<ounit(Mf,GeV)<<endl;*/ break;
    case 5 : Mf = interfaceBMass; /*cout<<"b"<<ounit(Mf,GeV)<<endl;*/ break;
  }
  if (Mf==0*GeV) 
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbar -- zero fermion mass."
		      << Exception::runerror; 
  double im=gw*Mf/2/MW;
  Complex c = Complex(0.,im);
  //cout<<"c: "<<c<<endl; 
  //cout<<"largeN wird berechnet werden";
  //cout<<"plusproduct"<<plusProduct(qbar,q)<<endl;
  //cout<<"minusproduct"<<minusProduct(qbar,q)<<endl;
  if (hel[qbar]==-hel[q]&& hel[q]==1){
    largeN = 0;
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==-hel[q]&& hel[q]==-1){
    largeN = 0;
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==hel[q] && hel[q]==1){
    largeN = c*(plusProduct(qbar,q));
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==hel[q] && hel[q]==-1){
    largeN = c*(minusProduct(qbar,q));
    //cout<<"largeN minus hat geklappt"<<largeN;
    return(largeN);
  }
  assert(false);
  return 0;
}

Complex MatchboxAmplitudehbbbar::evaluateOneLoop(size_t, const vector<int>& hel) {
  Complex res = 0;
  int q=0;
  int qbar=0;
  cPDVector x=amplitudePartonData();   
  
  
  for (;q<3;++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<3;++qbar){if (x[qbar]->id() ==-x[q]->id()) break;} 
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  
  long id=x[q]->id();
  Energy Mf = 0*GeV;
  switch(id){
    case 1 : Mf = interfaceDMass; break;
    case 2 : Mf = interfaceUMass; break;
    case 3 : Mf = interfaceSMass; break;
    case 4 : Mf = interfaceCMass; break;
    case 5 : Mf = interfaceBMass; break;
  }
  if (Mf==0*GeV) 
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbar -- zero fermion mass."
		      << Exception::runerror; 
  
  double loop = SM().alphaS()*CF/2/Constants::pi ; //one-loop-Factor
  double bornim = gw*Mf/2/MW; //constant factor from born
  double real = loop*bornim;
  Complex c = Complex(real,0.);
  
  if (hel[qbar]==-hel[q]&& hel[q]==1){
    res = 0;
    return(res);
  }
  if (hel[qbar]==-hel[q]&& hel[q]==-1){
    res = 0;
    return(res);
  }
  if (hel[qbar]==hel[q] && hel[q]==1){
    res = c*(plusProduct(qbar,q));
    return(res);
  }
  if (hel[qbar]==hel[q] && hel[q]==-1){
    res = c*(minusProduct(qbar,q));
    return(res);
  }
  assert(false);
  return 0;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehbbbar::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceUMass,GeV)<< ounit(interfaceDMass,GeV)<< ounit(interfaceSMass,GeV)<< ounit(interfaceCMass,GeV)<< ounit(interfaceBMass,GeV);
}
void MatchboxAmplitudehbbbar::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceUMass,GeV)>> iunit(interfaceDMass,GeV)>> iunit(interfaceSMass,GeV)>> iunit(interfaceCMass,GeV)>> iunit(interfaceBMass,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehbbbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehbbbar("Herwig::MatchboxAmplitudehbbbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehbbbar::Init() {

  static ClassDocumentation<MatchboxAmplitudehbbbar> documentation
    ("MatchboxAmplitudehbbbar");
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceUMass
    ("interfaceUMass",
     "The up quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceUMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceDMass
    ("interfaceDMass",
     "The down quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceDMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceSMass
    ("interfaceSMass",
     "The strange quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceSMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceCMass
    ("interfaceCMass",
     "The charm quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceCMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceBMass
    ("interfaceBMass",
     "The bottom quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceBMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  

}

