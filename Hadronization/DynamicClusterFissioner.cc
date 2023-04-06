// -*- C++ -*-
//
// ClusterFissioner.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// Thisk is the implementation of the non-inlined, non-templated member
// functions of the ClusterFissioner class.
//

#include "DynamicClusterFissioner.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<DynamicClusterFissioner,ClusterFissioner>
describeDynamicClusterFissioner("Herwig::DynamicClusterFissioner","Herwig.so");

DynamicClusterFissioner::DynamicClusterFissioner() :
  _Qqtilde(4*GeV),
  _Qg2tilde(4*GeV),
  _AngOrdFission(1),
  _restrictGluon(1)
{}

IBPtr DynamicClusterFissioner::clone() const {
  return new_ptr(*this);
}

IBPtr DynamicClusterFissioner::fullclone() const {
  return new_ptr(*this);
}

void DynamicClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << _AngOrdFission << _restrictGluon;;
}

void DynamicClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> _AngOrdFission >> _restrictGluon;;
}

void DynamicClusterFissioner::Init() {

    static Parameter<DynamicClusterFissioner,Energy> interfaceQqtilde
    ("Qqtilde",
     "Scale of the non-pert. quark splitting in the cluster fission",
     &DynamicClusterFissioner::_Qqtilde, GeV, 4.0*GeV, 1.5*GeV, 10000.0*GeV,false,false,false); //Set upper/lower bound and default  later

    static Parameter<DynamicClusterFissioner,Energy> interfaceQg2tilde
    ("Qg2tilde",
     "Upper scale for the Sudakov for the non-pert. gluon splitting in the cluster fission",
     &DynamicClusterFissioner::_Qg2tilde, GeV, 4.0*GeV, 0.5*GeV, 10000.0*GeV,false,false,false); //Set upper/lower bound and default later

  static Switch<DynamicClusterFissioner,int> interfaceAngOrdFission
    ("AngularOrderedFission",
     "Option to switch on angular ordering in the dynamic fission model",
     &DynamicClusterFissioner::_AngOrdFission,0,false,false);
  static SwitchOption interfaceAngOrdFissionNo
    (interfaceAngOrdFission,
     "No",
     "No (Qqtilde and Qg2tilde independent)",
     0);
  static SwitchOption interfaceAngOrdFissionYes
    (interfaceAngOrdFission,
     "Yes",
     "Qg2tilde=(1-z)*Qqtilde",
     1);


    static Switch<DynamicClusterFissioner,int> interfacerestrictGluon
    ("restrictGluon",
     "Option to restrict the gluon virtuality in the dynamic fission in order to make sure that the splitting works out kinematically",
     &DynamicClusterFissioner::_restrictGluon,0,false,false);
  static SwitchOption interfacerestrictGluonNo
    (interfacerestrictGluon,
     "No",
     "Gluon can have virtuality up to Qg2tilde. Might be necessary to start fission again with new z and qtilde for the quark",
     0);
  static SwitchOption interfacerestrictGluonYes
    (interfacerestrictGluon,
     "Yes",
     "Gluon virtuality is restricted to the range that the kinematics allow once that z and qitlde are fixed for the quark",
     1);

}




InvEnergy DynamicClusterFissioner::Pqzproposal(double z,Energy qtilde, Energy Qqtilde, Energy mq) const{
  return gluonMassGenerator()->ClusterAlphaS(0.*GeV2)*2.*(1.-sqr(mq/Qqtilde))/(qtilde*(1.-z)); 
}

InvEnergy DynamicClusterFissioner::Pqz(double z, Energy qtilde, Energy mq) const {
  if ( z*qtilde > mq){
    return gluonMassGenerator()->ClusterAlphaS(sqr(z*(1.0-z)*qtilde))*( (1.+z*z)/(1.-z) - ((2.*sqr(mq))/( z*(1.-z)*sqr(qtilde) )))/qtilde;}
  else {
    return 0./(1.*GeV);
  }
}


void DynamicClusterFissioner::dynamicFission(tPPtr & ptrP, tPPtr & ptrPbar, PPtr & ptrQ, PPtr & ptrQbar) const {

  
  LorentzMomentum P1, P2;
  bool fromcolor;

  Energy mq2, mg, MP1, Qqtilde, Qgtilde, qtilde, mq;

  double z;

  //minimal quark mass
  Energy m0, mu, ms, md;
  mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
  md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
  ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();
  m0=md;
  if(mu<m0){m0=mu;}
  if(ms<m0){m0=ms;}
  
 
  bool repeat2=true;


  LorentzMomentum Pcl(ptrP->momentum()+ptrPbar->momentum());
  Energy Mcl = Pcl.m();
  //boost to rest frame of the cluster
  Boost bv = Pcl.boostVector();


  //upper scale of the quark splitting
  Qqtilde = _Qqtilde;


 

while(repeat2){

    //choose from which constituent we radiate the gluon
      if ( UseRandom::rnd(0.0, 1.0)<0.5 ){
        P1 = ptrP->momentum();
        P2 = ptrPbar->momentum();
        mq = ptrP->data().constituentMass();
        mq2 = ptrPbar->data().constituentMass();
        fromcolor=true;
      }
      else{
        P2 = ptrP->momentum();
        P1 = ptrPbar->momentum();
        mq = ptrPbar->data().constituentMass();
        mq2 = ptrP->data().constituentMass();
        fromcolor=false;
      }

    

      //and check if it is kinematically possible for that constituent
      if (!canSplitMinimally(Mcl,mq,mq2,m0)){
        swap(P1,P2);
        swap(mq,mq2);
        fromcolor=(!fromcolor);
      }


   P1.boost(-bv);
  P2.boost(-bv);
  
  


  //get z and qtilde
  double zmin=mq/Qqtilde;
  double zmax=1-(4*sqr(m0)/( sqr(Mcl-mq2)-sqr(mq)  ));
  //if ((_AngOrdFission==1) && (1-(4*m0/Qqtilde)<zmax)){
 //   zmax=1-(4*m0/Qqtilde);
//  }
   

  

  double rzmin=-log(1-zmin);
  double rzmax=-log(1-zmax);
  
  bool repeat = true;
  while (repeat){
    //generate z and qtilde from a proposal distribution and check that at least with minimal quark mass the kinematics can work out
    bool repeat3 = true;
    double ymin, ymax;
    while (repeat3){
      z = 1.0 - exp(-UseRandom::rnd(rzmin,rzmax));

      ymin =mq/(z*Qqtilde);
      if ((_AngOrdFission==1)&&(4*m0/((1.-z)*Qqtilde)>ymin)){
        ymin=4*m0/((1.-z)*Qqtilde);
      }
      ymax =sqrt( ( sqr(Mcl-mq2)  -sqr(mq) -(4*sqr(m0)/(1-z)) )/( z*(1-z) ) )/Qqtilde;
      if (ymax >1.){
        ymax=1.;
      }
      if (ymin<ymax){
        repeat3=false;
      }
    }
  
      double rymin=-log(ymax);
      double rymax=-log(ymin);

      qtilde=Qqtilde*exp(-UseRandom::rnd(rymin,rymax));
     
    //compare with actual distrubition
    if (Pqzproposal(z,qtilde,Qqtilde,mq)*UseRandom::rnd(0.0,1.0)<Pqz(z,qtilde,mq)){
      repeat = false;
    }
  }
  
  

  //scale for the Sudakov in the gluon splitting
  if (_AngOrdFission==1){
    Qgtilde=(1.-z)*qtilde;
  }
  else{
    Qgtilde=_Qg2tilde;
  }
  

  //gluon virtuality
  if (_restrictGluon==1){
    //restirct the gluon virtuality to the range where kinematic reconstruction is possible, i.e. Mcl>MP1+mq2 is always true
    Energy mgmax=sqrt((1.-z)*(sqr(Mcl-mq2) - sqr(mq) -z*(1.-z)*sqr(qtilde)));
    mg=gluonMassGenerator()->generate(Qgtilde,mgmax);
  }
  else{
    //no additional restriction on gluon virtuality. Mcl>MP1+mq2 might fail and new z and qtilde for the quark are generated
    mg=gluonMassGenerator()->generate(Qgtilde);
  }

  

  



  //kinematic reconstruction / get momenta P1prime and P2prime before radiating
  //invariant mass of quark before radiating
   MP1 = sqrt(sqr(mq)+z*(1.0-z)*sqr(qtilde)+(sqr(mg)/(1.-z)));

 

  if(Mcl>MP1+mq2){
    repeat2 = false;
  }

 }//end while(repeat2)

  
 

  Energy absP = sqrt(sqr(sqr(Mcl))+sqr(sqr(MP1))+sqr(sqr(mq2))-2.0*(sqr(Mcl*MP1)+sqr(Mcl*mq2)+sqr(MP1*mq2)))/(2.0*Mcl);

  
 
  Lorentz5Momentum P1prime(absP*P1.vect().unit(),sqrt(sqr(absP)+sqr(MP1)),MP1);
  Lorentz5Momentum P2prime(absP*P2.vect().unit(),sqrt(sqr(absP)+sqr(mq2)),mq2);
 




  //splitting
  //construct the backwards light-like vector
  LorentzVector<double> nbar(-P1prime.vect().unit(),1.0);

  //construct the transverse momentum with random azimuthal angle
  //construct an orthonormal basis of transverse 4-vectors
  ThreeVector<double> Three_qperp1 = nbar.vect().orthogonal().unit();
  ThreeVector<double> Three_qperp2 = nbar.vect().unit().cross(Three_qperp1);
  LorentzVector<double> qperp1(Three_qperp1,0.0); 
  LorentzVector<double> qperp2(Three_qperp2,0.0);

 

  using Constants::pi;
  double phi = UseRandom::rnd( -pi , pi );
  Energy gperpabs=sqrt(sqr(1.-z)*(sqr(z*qtilde)-sqr(mq)) );

  LorentzMomentum gperp = gperpabs*cos(phi)*qperp1+gperpabs*sin(phi)*qperp2;

  //quark and gluon momenta
  Lorentz5Momentum k(z*P1prime + ((sqr(mq)-sqr(z*MP1)+sqr(gperpabs))/(2.0*z*P1prime.dot(nbar)))*nbar - gperp,mq);
  Lorentz5Momentum g((1.-z)*P1prime + ((sqr(mg)-sqr((1.-z)*MP1)+sqr(gperpabs))/(2.0*(1.-z)*P1prime.dot(nbar)))*nbar +gperp,mg);

  
 
  
  //create particles and give them to PartonSplitter
  

  PPtr ptrPprime = PPtr();
  PPtr ptrG = PPtr();
  ptrG = getParticle(ThePEG::ParticleID::g);
  ptrG->set5Momentum(g);
  if (fromcolor){
    ptrP->set5Momentum(k);
    ptrPbar->set5Momentum(P2prime);
    ptrPprime = getParticle(ptrP->id());
    ptrPprime->set5Momentum(P1prime);
    partonSplitter()->splitTimeLikeGluon(ptrG,ptrPprime,ptrPbar,ptrQ,ptrQbar,Qgtilde,nbar,fromcolor);
  }
  else{
    ptrP->set5Momentum(P2prime);
    ptrPbar->set5Momentum(k);
    ptrPprime = getParticle(ptrPbar->id());
    ptrPprime->set5Momentum(P1prime);
    partonSplitter()->splitTimeLikeGluon(ptrG,ptrP,ptrPprime,ptrQ,ptrQbar,Qgtilde,nbar,fromcolor);
  }




  //boost back to lab frame
  //(this must be possible in a nicer way)
  Lorentz5Momentum temp;
  temp = ptrP->momentum();
  temp.boost(bv);
  ptrP->set5Momentum(temp);

  temp = ptrPbar->momentum();
  temp.boost(bv);
  ptrPbar->set5Momentum(temp);

  temp = ptrQ->momentum();
  temp.boost(bv);
  ptrQ->set5Momentum(temp);

  temp = ptrQbar->momentum();
  temp.boost(bv);
  ptrQbar->set5Momentum(temp);


}


bool DynamicClusterFissioner::canSplitMinimally(Energy Mcl, Energy m1, Energy m2, Energy m0) const {

  if (Mcl<m1+m2+2*m0){ return false; }

  Energy Qtilde=_Qqtilde;


  if (Qtilde>m1){

  if (_AngOrdFission==0){
    if (Qtilde>2*m0+m1){
      return Mcl>m2+m1+2*m0;
    }
    else{
      return Mcl>m2+sqrt(Qtilde*(m1- (4*sqr(m0)/(Qtilde-m1))  ));
    }
  }
  else{
    if (Qtilde>4*m0+m1){
      return Mcl>m2+sqrt((m0+m1)*(4*m0+m1));
    }
    else{
      return Mcl<m2+sqrt( sqr(m1) + (4*sqr(m0)*(4*m1+Qtilde))/(Qtilde-m1)  );
    }
  }

  }
  else{
    return false;
  }


}

bool DynamicClusterFissioner::canSplitMinimally(tcClusterPtr clu, Energy minmass) const {
  Energy Mcl=clu->mass();
  Energy m1 = clu->particle(0)->data().constituentMass();
  Energy m2 = clu->particle(1)->data().constituentMass();
  return ((canSplitMinimally(Mcl,m1,m2,minmass))||(canSplitMinimally(Mcl,m2,m1,minmass)));
}
