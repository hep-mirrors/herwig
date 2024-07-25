// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GluonMassGenerator class.
//

#include "DynamicGluonMassGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ClusterHadronizationHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <ThePEG/Interface/Parameter.h>

using namespace Herwig;

DynamicGluonMassGenerator::DynamicGluonMassGenerator() {}

DynamicGluonMassGenerator::~DynamicGluonMassGenerator() {}

IBPtr DynamicGluonMassGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr DynamicGluonMassGenerator::fullclone() const {
  return new_ptr(*this);
}



InvEnergy DynamicGluonMassGenerator::PmgProposal(Energy, Energy mq) const {
    return (ClusterAlphaS(4*mq*mq)*0.32/mq);
}


InvEnergy DynamicGluonMassGenerator::PmgProposal(Energy mg) const {
  Energy mu, ms, md;
  mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
  md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
  ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();
  return PmgProposal(mg, mu)+PmgProposal(mg, md)+PmgProposal(mg, ms);
}



InvEnergy DynamicGluonMassGenerator::Pmg(Energy mg, Energy mq, Energy Qtilde) const {
  if((2*mq<mg) && (mg <= sqrt(mq*Qtilde))){
    return (ClusterAlphaS(mg*mg)/mg)*sqrt(1-pow(2*mq/mg,2))*(1+(2*pow(mq/mg,2)));
  }
  else if ((sqrt(mq*Qtilde)<mg) &&(mg<=Qtilde/2.0)){
    return (ClusterAlphaS(mg*mg)/mg)*sqrt(1-pow(2*mg/Qtilde,2))*(1+(3*pow(mq/mg,2))-pow(mg/Qtilde,2) );
  }
  else{
    return 0.0*InvGeV;
  }
}

InvEnergy DynamicGluonMassGenerator::Pmg(Energy mg,Energy Qtilde) const {
  Energy mu, ms, md;
  mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
  md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
  ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();
  return Pmg(mg,mu,Qtilde)+Pmg(mg,md,Qtilde)+Pmg(mg,ms,Qtilde);
}



Energy DynamicGluonMassGenerator::generateProposal(Energy mgmax) const {
  Energy m0, mu, ms, md;
  mu=getParticleData(ThePEG::ParticleID::u)->constituentMass();
  md=getParticleData(ThePEG::ParticleID::d)->constituentMass();
  ms=getParticleData(ThePEG::ParticleID::s)->constituentMass();
  m0=md;
  if(mu<m0){m0=mu;}
  if(ms<m0){m0=ms;}
  return UseRandom::rnd(2.*m0,mgmax);
}

Energy DynamicGluonMassGenerator::generate(Energy Qtilde, Energy mgmax) const {
  Energy mg;
  double r;
  bool repeat;
  Energy max=mgmax;
  if(0.5*Qtilde<max){max=0.5*Qtilde;};
    repeat=true;
    while(repeat){
        mg=generateProposal(max);
        r=UseRandom::rnd(0.0, 1.0);
        if(r*PmgProposal(mg)<Pmg(mg,Qtilde)){
          repeat=false;
        }
    }
    return mg;
}


Energy DynamicGluonMassGenerator::generate(Energy Qtilde) const {
  return generate(Qtilde,Qtilde);
}


Energy DynamicGluonMassGenerator::generate() const {
  return generate(Qgtilde());
}



// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



void DynamicGluonMassGenerator::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os  << ounit(_Qgtilde,GeV) << ounit(_clusteralphasfreeze,GeV);
}

void DynamicGluonMassGenerator::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> iunit(_Qgtilde,GeV) >> iunit(_clusteralphasfreeze,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DynamicGluonMassGenerator,GluonMassGenerator>
  describeHerwigDynamicGluonMassGenerator("Herwig::DynamicGluonMassGenerator", "");

void DynamicGluonMassGenerator::Init() {

  static ClassDocumentation<DynamicGluonMassGenerator> documentation
    ("There is no documentation for the DynamicGluonMassGenerator class");


  static Parameter<DynamicGluonMassGenerator,Energy> interfaceQgtilde
  ("Qgtilde",
   "Upper scale for the Sudakov for the non-pert. gluon splitting",
    &DynamicGluonMassGenerator::_Qgtilde, GeV, 4.0*GeV, 0.5*GeV, 10000.0*GeV,false,false,false);

  static Parameter<DynamicGluonMassGenerator,Energy> interfaceClusterAlphaSFreeze
  ("ClusterAlphaSFreeze",
   "Freeze-out scale for the non-pert. coupling in gluon splitting and cluster fission",
   &DynamicGluonMassGenerator::_clusteralphasfreeze, GeV, 1.0*GeV, 0.1*GeV, 10.0*GeV,false,false,false);

}

