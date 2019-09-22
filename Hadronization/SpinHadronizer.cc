// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinHadronizer class.
//

#include "SpinHadronizer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"
# include "Herwig/Utilities/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Cluster.h"

using namespace Herwig;

void SpinHadronizer::
handle(EventHandler &, const tPVector & tagged,const Hint & ) {
  for(const tPPtr & hadron : tagged) {
    // mesons
    if(MesonMatcher::Check(hadron->data())) {
      continue;
    }
    // baryons
    else if(BaryonMatcher::Check(hadron->data())) {
      baryonSpin(hadron);
    }
    else 
      continue;
  }
}

void SpinHadronizer::baryonSpin(tPPtr baryon) {
  // check only one parent
  if(baryon->parents().size()!=1) return;
  tPPtr parent = baryon->parents()[0];
  // and its a cluster
  if(parent->id()!=ParticleID::Cluster) return;
  tClusterPtr cluster = dynamic_ptr_cast<tClusterPtr>(parent);
  int prim_quark = (abs(baryon->id())/1000)%10;
  int sign_quark = baryon->id()>0 ? prim_quark : -prim_quark;
  // only strange, charm and bottom for the moment
  if(prim_quark<3) return;
  tPPtr quark;
  for(unsigned int ix=0;ix<cluster->numComponents();++ix) {
    if(cluster->particle(ix)->id()==sign_quark) {
      quark = cluster->particle(ix);
    }
  }
  if(!quark) return;
  if(!quark->spinInfo()) return;
  tcFermionSpinPtr sp(dynamic_ptr_cast<tcFermionSpinPtr>(quark->spinInfo()));
  // decay it
  sp->decay();
  // create the spin info
  if(baryon->dataPtr()->iSpin()==PDT::Spin1Half) {
    vector<SpinorWaveFunction> waves;
    RhoDMatrix rho;
    SpinorWaveFunction::calculateWaveFunctions(waves,rho,baryon,outgoing);
    SpinorWaveFunction::constructSpinInfo(waves,baryon,outgoing,true);
  }
  else if(baryon->dataPtr()->iSpin()==PDT::Spin3Half) {
    vector<RSSpinorWaveFunction> waves;
    RhoDMatrix rho;
    RSSpinorWaveFunction::calculateWaveFunctions(waves,rho,baryon,outgoing);
    RSSpinorWaveFunction::constructSpinInfo(waves,baryon,outgoing,true);
  }
  // can't handle spin 5/2 > 3/2
  else {
    return;
  }
  // extract the polarization of the quark
  double pol = 2.*sp->rhoMatrix()(1,1).real()-1.;
  // the different options for different spin types
  vector<double> polB(baryon->dataPtr()->iSpin(),1./double(baryon->dataPtr()->iSpin()));
  const int mult = prim_quark*1000;
  int bid = abs(baryon->id());
  // lambda and Xi spin 1/2 (spin0 diquark)
  if(bid== mult+122|| bid== mult+132|| bid== mult+232) {
    baryon->spinInfo()->rhoMatrix()(0,0) = 0.5*(1.-pol);
    baryon->spinInfo()->rhoMatrix()(1,1) = 0.5*(1.+pol);
  }
  // sigma_b, xi' and omega_b spin 1/2 (spin1 diquark)
  else if(bid== mult+112|| bid== mult+212|| bid== mult+222||
	  bid== mult+312|| bid== mult+322|| bid== mult+332) {
    baryon->spinInfo()->rhoMatrix()(0,0) = 0.5*(1.-pol) +pol*omegaHalf_;
    baryon->spinInfo()->rhoMatrix()(1,1) = 0.5*(1.+pol) -pol*omegaHalf_;
  }
  // sigma*, xi* and omegab* spin 3/2 (spin1 diquark)
  else if(bid== mult+114|| bid== mult+214|| bid== mult+224|| bid== mult+334) { 
    baryon->spinInfo()->rhoMatrix()(0,0) = 0.375*(1.-pol)*omegaHalf_;
    baryon->spinInfo()->rhoMatrix()(1,1) = 0.5*(1.-pol)-omegaHalf_/6.*(3.-5.*pol);
    baryon->spinInfo()->rhoMatrix()(2,2) = 0.5*(1.+pol)-omegaHalf_/6.*(3.+5.*pol);
    baryon->spinInfo()->rhoMatrix()(3,3) = 0.375*(1.+pol)*omegaHalf_;
  }
  else
    return;
  

  
  // generator()->log() << "Baryon: " << *baryon << "\n";
  // generator()->log() << "Parent: " << *cluster << "\n";
  // generator()->log() << "Quark: " << *quark << "\n";
  // generator()->log() << "Rho\n" << sp->rhoMatrix() << "\n";
  // generator()->log() << "testing is decayed " << sp->decayed() <<" \n";
  // generator()->log() << baryon->spinInfo()->rhoMatrix() << "\n";
}

IBPtr SpinHadronizer::clone() const {
  return new_ptr(*this);
}

IBPtr SpinHadronizer::fullclone() const {
  return new_ptr(*this);
}

void SpinHadronizer::persistentOutput(PersistentOStream & os) const {
  os << omegaHalf_;
}

void SpinHadronizer::persistentInput(PersistentIStream & is, int) {
  is >> omegaHalf_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SpinHadronizer,StepHandler>
  describeHerwigSpinHadronizer("Herwig::SpinHadronizer", "Herwig.so");

void SpinHadronizer::Init() {

  static ClassDocumentation<SpinHadronizer> documentation
    ("The SpinHadronizer class implements a simple mode for"
     " the transfer of spin from quarks to hadrons");

  static Parameter<SpinHadronizer,double> interfaceOmegaHalf
    ("OmegaHalf",
     "The omega_1/2 Falk-Psekin parameter",
     &SpinHadronizer::omegaHalf_, 2./3., 0.0, 1.0,
     false, false, Interface::limited);

}
