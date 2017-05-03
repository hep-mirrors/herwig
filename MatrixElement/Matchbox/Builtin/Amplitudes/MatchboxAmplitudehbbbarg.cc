// -*- C++ -*-
//
// MatchboxAmplitudehbbbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehbbbarg class.
//

#include "MatchboxAmplitudehbbbarg.h"
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

MatchboxAmplitudehbbbarg::MatchboxAmplitudehbbbarg() {}

MatchboxAmplitudehbbbarg::~MatchboxAmplitudehbbbarg() {}

IBPtr MatchboxAmplitudehbbbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehbbbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehbbbarg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);;
}

void MatchboxAmplitudehbbbarg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehbbbarg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  /*for (PDVector::iterator x=xproc.begin(); x!=xproc.end();  ++x){
    cout<<"xproc: "<<(**x).id()<<endl;
  }
  cout<<endl;*/
  PDVector::iterator q=xproc.begin();
  for (; q!=xproc.end(); ++q){
      if ((**q).id()==1 || 
          (**q).id()==2 ||
          (**q).id()==3 ||
          (**q).id()==4 ||
          (**q).id()==5) 
          break;
      //cout<<"kein Quark gefunden";
  }
  if(q==xproc.end()) return false;
  int qid = (**q).id(); //cout<<"qid ="<<qid<<endl;
  xproc.erase(q);
  PDVector::iterator qb=xproc.begin();
  for (; qb!=xproc.end(); ++qb){
      if ((**qb).id()==-qid /*|| (**qb).id()==qid*/ ){/*cout<<"qbid = "<<(**qb).id()<<endl; cout<<"kein antiquark gefunden!"<<endl;*/ break;}
  }
  if(qb==xproc.end()){/*cout<<"kein antiquark gefunden"<<endl;*/ return false;}
  xproc.erase(qb);
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21){xproc.erase(g); break; }
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehbbbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  /*setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));*/

  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehbbbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  unsigned int g=0;
  cPDVector x=amplitudePartonData(); 
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;}
  for (;g<amplitudePartonData().size();++g){if (x[g]->id() == 21) break;}
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  //double alphaS = SM().alphaS();
  //double v= 2*MW/gw/sqrt(mu2()) ; //Auf mu normierter VakuumErwartungswert
  //double c = gs*alphaS/3/Constants::pi/v;
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
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbarg -- zero fermion mass."
		      << Exception::runerror; 
  Complex c2= Complex(0.,gw*gs*Mf/MW/sqrt(2));
  
  //qghq Prozesse mit Quark als "Zwischenteilchen" aus qghq.nb
  if(hel[qbar]==+1 && hel[q]==+1 && hel[g]==-1){
    largeN = -c2*minusProduct(q,qbar)*minusProduct(q,qbar)/minusProduct(q,g)/minusProduct(qbar,g);
    return(largeN);
  }
  if(hel[qbar]==+1 && hel[q]==+1 && hel[g]==+1){
    largeN = c2*(minusProduct(qbar,g)/plusProduct(q,g)+
                 minusProduct(q,g)/plusProduct(qbar,g)-
                 invariant(qbar,q)/plusProduct(q,g)/plusProduct(qbar,g));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==-1 && hel[g]==+1){
    largeN = c2*plusProduct(q,qbar)*plusProduct(q,qbar)/plusProduct(q,g)/plusProduct(qbar,g);
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==-1 && hel[g]==-1){
    largeN = -c2*(plusProduct(qbar,g)/minusProduct(q,g)+
                  plusProduct(q,g)/minusProduct(qbar,g)-
                  invariant(qbar,q)/minusProduct(q,g)/minusProduct(qbar,g));
    return(largeN);
  }
  if(hel[qbar]==-hel[q]){
    largeN = 0;
    return(largeN);
  }
  
//  Viel zu hohe Ordnnung in alphaS!!! Vielleicht spaeter per Interface dazuschaltbar??  
//  //mpm verdrehte vorzeichen wegen simons spinordef q=p qbar=q; g=r gluonhelizitaeten bleiben bleiben
//  if(hel[qbar]==+1 && hel[q]==-1 && hel[g]==-1){
//      largeN = c*minusProduct(q,qbar)*plusProduct(q,g)*plusProduct(q,g)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //mpp
//  if(hel[qbar]==1 && hel[q]==-1 && hel[g]==1){
//      largeN = c*minusProduct(qbar,g)*minusProduct(qbar,g)*plusProduct(q,qbar)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //pmm
//  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==-1){
//      largeN = -c*minusProduct(q,qbar)*plusProduct(qbar,g)*plusProduct(qbar,g)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //pmp
//  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==1){
//      largeN = -c*minusProduct(q,g)*minusProduct(q,g)*plusProduct(q,qbar)/(invariant(q,qbar));
//      return(largeN);
//  }
  assert(false);

  return 0;
}

/*Complex MatchboxAmplitudehbbbarg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehbbbarg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceUMass,GeV)<< ounit(interfaceDMass,GeV)<< ounit(interfaceSMass,GeV)<< ounit(interfaceCMass,GeV)<< ounit(interfaceBMass,GeV);
}
void MatchboxAmplitudehbbbarg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceUMass,GeV)>> iunit(interfaceDMass,GeV)>> iunit(interfaceSMass,GeV)>> iunit(interfaceCMass,GeV)>> iunit(interfaceBMass,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehbbbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehbbbarg("Herwig::MatchboxAmplitudehbbbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehbbbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudehbbbarg> documentation
    ("MatchboxAmplitudehbbbarg");
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceUMass
    ("interfaceUMass",
     "The up quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceUMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceDMass
    ("interfaceDMass",
     "The down quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceDMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceSMass
    ("interfaceSMass",
     "The strange quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceSMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceCMass
    ("interfaceCMass",
     "The charm quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceCMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceBMass
    ("interfaceBMass",
     "The bottom quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceBMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  

}

