// -*- C++ -*-
//
// PartonSplitter.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonSplitter class.
//

#include "DynamicPartonSplitter.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/EventRecord/Step.h>
#include "ThePEG/Interface/Parameter.h"
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include "ClusterHadronizationHandler.h"
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/PDT/PDT.h>


using namespace Herwig;



void DynamicPartonSplitter::persistentOutput(PersistentOStream & os) const {
  os << _findProgenitor << _restrictZ << _dynamicGluonMassGenerator;
}

void DynamicPartonSplitter::persistentInput(PersistentIStream & is, int) {
  is >> _findProgenitor >> _restrictZ >> _dynamicGluonMassGenerator;
}

DescribeClass<DynamicPartonSplitter,PartonSplitter>
describeDynamicPartonSplitter("Herwig::DynamicPartonSplitter","Herwig.so");


void DynamicPartonSplitter::Init() {

  static ClassDocumentation<DynamicPartonSplitter> documentation
    ("This class is reponsible of the nonperturbative splitting of partons");


  static Switch<DynamicPartonSplitter,int> interfaceRestrictZ
    ("RestrictZ",
     "Option to restrict z in gluon splitting to z > 0.5",
     &DynamicPartonSplitter::_restrictZ,0,false,false);
  static SwitchOption interfaceRestrictZNo
    (interfaceRestrictZ,
     "No",
     "No restriction",
     0);
  static SwitchOption interfaceRestrictZYes
    (interfaceRestrictZ,
     "Yes",
     "Restric z > 0.5",
     1);


  static Switch<DynamicPartonSplitter,int> interfaceFindProgenitor
    ("FindProgenitor",
     "Option for how to choose which color partner is the progenitor of the gluon",
     &DynamicPartonSplitter::_findProgenitor,0,false,false);
  static SwitchOption interfaceFindProgenitorAngle
    (interfaceFindProgenitor,
     "Angle",
     "choose color partner with maximum cos(theta)",
     0);
  static SwitchOption interfaceFindProgenitorMomentum
    (interfaceFindProgenitor,
     "Momentum",
     "choose color partner with maximum |p|*cos(theta)",
     1);
    static SwitchOption interfaceFindProgenitorpT
    (interfaceFindProgenitor,
     "pT",
     "choose color partner with maximum cos^2(theta) as backwards direction",
     2);
    static SwitchOption interfaceFindProgenitorpTforward
    (interfaceFindProgenitor,
     "pTforward",
     "choose color partner with maximum cos^2(theta), gluon always in forward direction",
     3);
    static SwitchOption interfaceFindProgenitorGluon
    (interfaceFindProgenitor,
     "Gluon",
     "Gluon as progenitor",
     4);
    static SwitchOption interfaceFindProgenitorIsotropic
    (interfaceFindProgenitor,
     "Isotropic",
     "Isotropic gluon decay",
     5);

    static Reference<DynamicPartonSplitter,DynamicGluonMassGenerator> interfaceDynamicGluonMassGenerator
    ("GluonMassGenerator",
     "Set a reference to a gluon mass generator.",
     &DynamicPartonSplitter::_dynamicGluonMassGenerator, false, false, true, true, false);
}


void DynamicPartonSplitter::drawNewFlavour(PPtr & ptrQ, PPtr & ptrQbar, Energy mg, Energy Qtilde){
  static const Energy md = getParticleData(ThePEG::ParticleID::d)->constituentMass();
  static const Energy mu = getParticleData(ThePEG::ParticleID::u)->constituentMass();
  static const Energy ms = getParticleData(ThePEG::ParticleID::s)->constituentMass();

  double prob_d=_dynamicGluonMassGenerator->Pmg(mg,md,Qtilde)*UnitRemoval::E;
  double prob_u=_dynamicGluonMassGenerator->Pmg(mg,mu,Qtilde)*UnitRemoval::E;
  double prob_s=_dynamicGluonMassGenerator->Pmg(mg,ms,Qtilde)*UnitRemoval::E;
 

  long idNew=0;
  int choice = UseRandom::rnd3(prob_u, prob_d, prob_s);
  switch(choice) {
    case 0: idNew = ThePEG::ParticleID::u; break;
    case 1: idNew = ThePEG::ParticleID::d; break;
    case 2: idNew = ThePEG::ParticleID::s; break;
  }
  ptrQ = getParticle(idNew);
  ptrQbar = getParticle(-idNew);
  assert (ptrQ);
  assert(ptrQbar);
  assert (ptrQ->dataPtr());
  assert(ptrQbar->dataPtr());

}


void DynamicPartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon, PPtr & ptrQ, PPtr & ptrQbar, tcPPtr ptrP, tcPPtr ptrPbar, Energy Qtilde, LorentzVector<double> nbar, bool ProgenitorHasColor){

  LorentzMomentum momentumG = ptrGluon->momentum();
  LorentzMomentum momentumP = ptrP->momentum();
  LorentzMomentum momentumPbar = ptrPbar->momentum();
  Energy mg = momentumG.m();

  drawNewFlavour(ptrQ,ptrQbar,mg,Qtilde);

  Energy mQ = ptrQ->data().constituentMass();
  LorentzMomentum momentumQ, momentumQbar;

  if(_findProgenitor == 5){

    Lorentz5Momentum momentum5Qtemp, momentum5Qbartemp;
    double cosThetaStar = UseRandom::rnd( -1.0 , 1.0 );
    using Constants::pi;
    double phiStar = UseRandom::rnd( -pi , pi );
    Kinematics::twoBodyDecay(momentumG, mQ,
         mQ, cosThetaStar, phiStar, momentum5Qtemp,
         momentum5Qbartemp );
    momentumQ = momentum5Qtemp;
    momentumQbar = momentum5Qbartemp;
   }
   else{
    //construct an orthonormal basis of transverse 4-vectors
    ThreeVector<double> Three_qperp1 = nbar.vect().orthogonal().unit();
    ThreeVector<double> Three_qperp2 = nbar.vect().unit().cross(Three_qperp1);
    LorentzVector<double> qperp1(Three_qperp1,0.0); 
    LorentzVector<double> qperp2(Three_qperp2,0.0);

    //get the z-value from the gluon splitting function
    bool repeat;
    double z,r,x;
    if (mg < sqrt(mQ*Qtilde)){x = pow(mQ/mg,2);}
      else {x = pow(mg/Qtilde,2);}
    repeat = true;
    while(repeat){
      r = UseRandom::rnd(0.0, 1.0);
      z = (r-0.5)*sqrt(1-4*x)+0.5;
      r = UseRandom::rnd(0.0, 1.0-2.0*x+2.0*pow(mQ/mg,2));
      if( r < (1-(2*z*(1-z))+2*pow(mQ/mg,2)) ){
        repeat = false;
      }
    }

    if (_restrictZ==1){//without this ProgenitorHasColor would not be necessary
      if (ProgenitorHasColor){
        if (z>0.5){z = 1.-z;} //anti-quark in direction of progenitor
      }
      else{
        if (z<0.5){z = 1.-z;} //quark in direction of progenitor
      }
    } 
 
    //magniuted of relative transverse momentum
    Energy qT = sqrt(sqr(mg)*z*(1.-z)-sqr(mQ));
    //azimuthal angle of relative tranverse momentum
    using Constants::pi;
    double phi = UseRandom::rnd( -pi , pi );

    //set the relative transverse momentum
    LorentzMomentum qperp = qT*cos(phi)*qperp1+qT*sin(phi)*qperp2;

    //calculate the quarks' momenta
    momentumQ = z*momentumG + ( ( (mg*mg)*(1-2.*z) -2*momentumG.dot(qperp)  )/( 2.*momentumG.dot(nbar) ) )*nbar + qperp;
    momentumQbar = momentumG - momentumQ;
 
  } //else _findProgenitor Isotropic

  //without this I would not need to give also P and Pbar
  if (((_findProgenitor == 4)||(_findProgenitor==5))&&(_restrictZ==1)){
    LorentzMomentum c1 = momentumQ+momentumP;
    LorentzMomentum c2 = momentumQbar+momentumPbar;
    LorentzMomentum c3 = momentumQ+momentumPbar;
    LorentzMomentum c4 = momentumQbar+momentumP;

    //linear
    Energy M12=c1.m()+c2.m();
    Energy M34=c3.m()+c4.m();
    
    //quadratic
    //Energy2 M12 = c1.m()*c1.m() + c2.m()*c2.m();
    //Energy2 M34 = c3.m()*c3.m() + c4.m()*c4.m();

    //product
    //Energy2 M12=c1.m()*c2.m();
    //Energy2 M34=c3.m()*c4.m();

    //maximum
    //Energy M12 = std::max(c1.m(),c2.m());
    //Energy M34 = std::max(c3.m(),c4.m());

    if (M12 < M34){
      swap(momentumQ,momentumQbar);
    }
  }

  Lorentz5Momentum momentum5Q(momentumQ,mQ);
  Lorentz5Momentum momentum5Qbar(momentumQbar,mQ);
  
  ptrQ    ->set5Momentum( momentum5Q    );
  ptrQbar ->set5Momentum( momentum5Qbar );

}



void DynamicPartonSplitter::splitTimeLikeGluon(tcPPtr ptrGluon, PPtr & ptrQ, PPtr & ptrQbar){
  LorentzVector<double> nbar;
  bool ProgenitorHasColor;
      
  tPPtr ColorPartner = ptrGluon->antiColourLine()->endParticle();
  tPPtr AntiColorPartner = ptrGluon->colourLine()->startParticle();
  Findnbar(ptrGluon,ColorPartner,AntiColorPartner,nbar,ProgenitorHasColor);
      
  splitTimeLikeGluon(ptrGluon,ptrQ,ptrQbar,ColorPartner,AntiColorPartner,_dynamicGluonMassGenerator->Qgtilde(),nbar,ProgenitorHasColor);
}



void DynamicPartonSplitter::Findnbar(tcPPtr ptrGluon, tcPPtr ptrP, tcPPtr ptrPbar, LorentzVector<double> & nbarfinal, bool & ProgenitorHasColor){

  LorentzMomentum momentumG = ptrGluon->momentum();
  LorentzMomentum momentumP = ptrP->momentum();
  LorentzMomentum momentumPbar = ptrPbar->momentum();

  LorentzMomentum nbar = momentumP;
  ProgenitorHasColor = true;
  
  if (_findProgenitor==0){
    if( momentumG.vect().unit().dot( momentumPbar.vect().unit() ) > momentumG.vect().unit().dot( momentumP.vect().unit() )){
      nbar = momentumPbar; //minimum angle theta
      ProgenitorHasColor = false; //has anticolor
    }
     nbar.setVect( (-1.)*nbar.vect() );
  }
  else if (_findProgenitor==1){
    if ( momentumG.vect().dot( momentumPbar.vect() ) > momentumG.vect().dot( momentumP.vect() ) ){ 
      nbar = momentumPbar; //maximum |p|*cos(theta)
      ProgenitorHasColor = false; //has anticolor
    }
     nbar.setVect( (-1.)*nbar.vect() );
  }
  else if ((_findProgenitor==2)||(_findProgenitor==3)){
    ProgenitorHasColor = false;
    if( pow(momentumG.vect().unit().dot( momentumPbar.vect().unit()),2) > pow(momentumG.vect().unit().dot( momentumP.vect().unit() ),2)){
      nbar = momentumPbar; //maximum cos^2(theta) --> minimum gluon transverse momentum
      ProgenitorHasColor = true;
    }
    if (_findProgenitor==3){
      if ( momentumG.vect().unit().dot( nbar.vect().unit() ) > 0.  ){
        nbar.setVect( (-1.)*nbar.vect() ); //gluon always in forward direction
        ProgenitorHasColor = !ProgenitorHasColor;
      }
    }
  }
  else if (_findProgenitor==4){
    nbar = momentumG;
    nbar.setVect( (-1.)*nbar.vect() ); 
  }
  
  nbarfinal.setT(1.0);
  nbarfinal.setVect(nbar.vect().unit());
}


