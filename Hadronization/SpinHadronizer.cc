// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinHadronizer class.
//

#include "SpinHadronizer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Cluster.h"

using namespace Herwig;

// choose only baryons
void SpinHadronizer::
handle(EventHandler &, const tPVector & tagged,const Hint & ) {
  for(const tPPtr & hadron : tagged) {
    // mesons
    if(MesonMatcher::Check(hadron->data())) {
      mesonSpin(hadron);
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
  unsigned int prim_quark = (abs(baryon->id())/1000)%10;
  int sign_quark = baryon->id()>0 ? prim_quark : -prim_quark;
  // only strange, charm and bottom for the moment
  if(prim_quark<minFlav_ || prim_quark>maxFlav_ ) return;
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
    SpinorWaveFunction::calculateWaveFunctions(waves,rho,baryon,baryon->id() >0 ? incoming : outgoing);
    SpinorWaveFunction::constructSpinInfo(waves,baryon,outgoing,true);
  }
  else if(baryon->dataPtr()->iSpin()==PDT::Spin3Half) {
    vector<RSSpinorWaveFunction> waves;
    RhoDMatrix rho;
    RSSpinorWaveFunction::calculateWaveFunctions(waves,rho,baryon,baryon->id() >0 ? incoming : outgoing);
    RSSpinorWaveFunction::constructSpinInfo(waves,baryon,outgoing,true);
  }
  // can't handle spin > 3/2
  else {
    return;
  }
  // extract the polarization of the quark
  double pol = 2.*sp->rhoMatrix()(1,1).real()-1.;
  if(sign_quark<0) {
    qPol_[prim_quark-3].first  += pol;
    qPol_[prim_quark-3].second += 1.;
  }
  else {
    qPol_[prim_quark].first  += pol;
    qPol_[prim_quark].second += 1.;
  }
  // the different options for different spin types
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
    baryon->spinInfo()->rhoMatrix()(1,1) = 0.5*(1.-pol)-omegaHalf_/8.*(3.-5.*pol);
    baryon->spinInfo()->rhoMatrix()(2,2) = 0.5*(1.+pol)-omegaHalf_/8.*(3.+5.*pol);
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

void SpinHadronizer::mesonSpin(tPPtr meson) {
  // check only one parent
  if(meson->parents().size()!=1) return;
  tPPtr parent = meson->parents()[0];
  // and its a cluster
  if(parent->id()!=ParticleID::Cluster) return;
  tClusterPtr cluster = dynamic_ptr_cast<tClusterPtr>(parent);
  unsigned int prim_quark = (abs(meson->id())/100)%10;
  // find the quark id
  int sign_quark = meson->id()>0 ? prim_quark : -prim_quark;
  // B and K mesons have antiquarks bbar,sbar in mesons and quarks, b,s, in antimesons
  if(prim_quark%2==1) sign_quark *= -1.;
  // only strange, charm and bottom for the moment
  if(prim_quark<minFlav_ || prim_quark>maxFlav_ ) return;
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
  if(meson->dataPtr()->iSpin()==PDT::Spin0) {
    RhoDMatrix rho;
    ScalarWaveFunction::calculateWaveFunctions(rho,meson,meson->id() > 0 ? incoming : outgoing);
    ScalarWaveFunction::constructSpinInfo(meson,outgoing,true);
  }
  else if(meson->dataPtr()->iSpin()==PDT::Spin1) {
    vector<VectorWaveFunction> waves;
    RhoDMatrix rho;
    VectorWaveFunction::calculateWaveFunctions(waves,rho,meson,meson->id() > 0 ? incoming : outgoing,false);
    VectorWaveFunction::constructSpinInfo(waves,meson,outgoing,true,false);
  }
  else if(meson->dataPtr()->iSpin()==PDT::Spin2) {
    vector<TensorWaveFunction> waves;
    RhoDMatrix rho;
    TensorWaveFunction::calculateWaveFunctions(waves,rho,meson,meson->id() > 0 ? incoming : outgoing,false);
    TensorWaveFunction::constructSpinInfo(waves,meson,outgoing,true,false);
  }
  else {
    return;
  }
  // extract the polarization of the quark
  double pol = 2.*sp->rhoMatrix()(1,1).real()-1.;
  if(sign_quark<0) {
    qPol_[prim_quark-3].first  += pol;
    qPol_[prim_quark-3].second += 1.;
  }
  else {
    qPol_[prim_quark].first  += pol;
    qPol_[prim_quark].second += 1.;
  }
  // the different options for different spin types
  int bid = abs(meson->id());
  // light-quark spin 1/2+ -> spin 0 heavy meson
  if(bid==10313 || bid==10323 || bid==10333 || bid==10413 || bid==10423
                || bid==10433 || bid==20413 || bid==20423 || bid==20433
                || bid==10513 || bid==10523 || bid==10533) {
    // Falk-Peskin "no-win" theorem for non-excited heavy mesons:
    // no polarization information would be find in the non-excited meson
    // for the excted mesons
    meson->spinInfo()->rhoMatrix()(0,0) = (1.-pol)/16. + (omegaThreeHalf_/16.)*(3.-5.*pol);
    meson->spinInfo()->rhoMatrix()(1,1) = 0.25*(1.-omegaThreeHalf_);
    meson->spinInfo()->rhoMatrix()(2,2) = (1.+pol)/16. + (omegaThreeHalf_/16.)*(3.+5.*pol);
  }
  // light-quark spin 3/2+ -> exited spin 2 meson
  else if(bid==315 || bid==325 || bid==335 || bid==415 || bid==425  || bid==435
                   || bid==515 || bid==525 || bid==535) {
    meson->spinInfo()->rhoMatrix()(0,0) = 0.25*(1.-pol)*omegaThreeHalf_;
    meson->spinInfo()->rhoMatrix()(1,1) = 0.1875*(1.-pol)-0.125*(1.-pol)*omegaThreeHalf_;
    meson->spinInfo()->rhoMatrix()(2,2) = 0.25*(1.-omegaThreeHalf_);
    meson->spinInfo()->rhoMatrix()(3,3) = 0.1875*(1.+pol)-0.125*(1.+pol)*omegaThreeHalf_;
    meson->spinInfo()->rhoMatrix()(4,4) = 0.25*(1.+pol)*omegaThreeHalf_;
  }
  else {
    return;
  }
}


IBPtr SpinHadronizer::clone() const {
  return new_ptr(*this);
}

IBPtr SpinHadronizer::fullclone() const {
  return new_ptr(*this);
}

void SpinHadronizer::persistentOutput(PersistentOStream & os) const {
  os << omegaHalf_ << omegaThreeHalf_;
}

void SpinHadronizer::persistentInput(PersistentIStream & is, int) {
  is >> omegaHalf_ >> omegaThreeHalf_;
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
     "The omega_1/2 Falk-Peskin parameter",
     &SpinHadronizer::omegaHalf_, 2./3., 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<SpinHadronizer,double> interfaceOmegaThreeHalf
    ("OmegaThreeHalf",
     "The omega_3/2 Falk-Peskin parameter",
     &SpinHadronizer::omegaThreeHalf_, 0.2, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<SpinHadronizer,unsigned int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of quark for which to transfer the polarization",
     &SpinHadronizer::minFlav_, 3, 3, 5,
     false, false, Interface::limited);

  static Parameter<SpinHadronizer,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of quark for which to transfer the polarization",
     &SpinHadronizer::maxFlav_, 5, 3, 5,
     false, false, Interface::limited);

  static Switch<SpinHadronizer,bool> interfaceDebug
    ("Debug",
     "Output info on polarizations each for debugging",
     &SpinHadronizer::debug_, false, false, false);
  static SwitchOption interfaceDebugYes
    (interfaceDebug,
     "Yes",
     "Debug",
     true);
  static SwitchOption interfaceDebugNo
    (interfaceDebug,
     "No",
     "No info",
     false);

}

void SpinHadronizer::doinit() {
  StepHandler::doinit();
  if(minFlav_>maxFlav_)
    throw InitException() << "The minimum flavour " << minFlav_
			  << "must be lower the than maximum flavour " << maxFlav_
			  << " in SpinHadronizer::doinit() "
			  << Exception::runerror;
}

void SpinHadronizer::dofinish() {
  StepHandler::dofinish();
  if(debug_) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(qPol_[ix].second!=0)
	generator()->log() << "Average polarization of " << getParticleData(long(3+ix))->PDGName() << " antiquarks "
			   << qPol_[ix].first/qPol_[ix].second << "\n";
      if(qPol_[ix+3].second!=0)
	generator()->log() << "Average polarization of " << getParticleData(long(3+ix))->PDGName()    << "     quarks "
			   << qPol_[ix+3].first/qPol_[ix+3].second << "\n";
    }
  }
}

void SpinHadronizer::doinitrun() {
  StepHandler::doinitrun();
}
