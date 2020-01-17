// -*- C++ -*-
//
// FFVCurrentDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVCurrentDecayer class.
//

#include "FFVCurrentDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

IBPtr FFVCurrentDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVCurrentDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVCurrentDecayer::doinit() {
  FFVPtr_ = dynamic_ptr_cast<FFVVertexPtr>(vertex());
  GeneralCurrentDecayer::doinit();
}


void FFVCurrentDecayer::rebind(const TranslationMap & trans)
  {
  FFVPtr_ = trans.translate(FFVPtr_);
  GeneralCurrentDecayer::rebind(trans);
}

IVector FFVCurrentDecayer::getReferences() {
  IVector ret = GeneralCurrentDecayer::getReferences();
  ret.push_back(FFVPtr_);
  return ret;
}

void FFVCurrentDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFVPtr_;
}

void FFVCurrentDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFVPtr_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFVCurrentDecayer,GeneralCurrentDecayer>
describeHerwigFFVCurrentDecayer("Herwig::FFVCurrentDecayer", "Herwig.so");

void FFVCurrentDecayer::Init() {

  static ClassDocumentation<FFVCurrentDecayer> documentation
    ("There is no documentation for the FFVCurrentDecayer class");

}

double FFVCurrentDecayer::me2(const int ichan, const Particle & inpart,
			      const ParticleVector & decay,
			      MEOption meopt) const {
  // get the particles for the hadronic curret
  Energy q;
  ParticleVector hadpart(decay.begin()+1,decay.end());
  // fermion types
  int itype[2];
  if(inpart.dataPtr()->CC())    itype[0] = inpart.id() > 0 ? 0 : 1;
  else                          itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                          itype[1] = 2;
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&inpart),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    // for the decaying particle
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(wave_,decay[0],outgoing,true);
    }
    weakCurrent()->current(mode(),ichan,q,hadpart,meopt);
    return 0.;
  }
  Energy2 scale(sqr(inpart.mass()));
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,decay[0],outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(wave_   ,decay[0],outgoing);
  // calculate the hadron current
  vector<LorentzPolarizationVectorE> 
    hadron(weakCurrent()->current(mode(),ichan,q,hadpart,meopt));
  // prefactor
  double pre = sqr(pow(inpart.mass()/q,int(hadpart.size()-2)));
  // work out the mapping for the hadron vector
  vector<unsigned int> constants(decay.size()+1),ihel(decay.size()+1);
  vector<PDT::Spin> ispin(decay.size());
  int itemp(1);
  unsigned int hhel,ix(decay.size());
  do {
    --ix;
    ispin[ix]=decay[ix]->data().iSpin();
    itemp*=ispin[ix];
    constants[ix]=itemp;
  }
  while(ix>0);
  constants[decay.size()]=1;
  constants[0]=constants[1];
  // compute the matrix element
  GeneralDecayMEPtr newME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,ispin)));
  VectorWaveFunction vWave;
  tcPDPtr vec= inpart.dataPtr()->iCharge()-decay[0]->dataPtr()->iCharge() > 0
    ? getParticleData(ParticleID::Wplus) : getParticleData(ParticleID::Wminus);
  Lorentz5Momentum vmom=inpart.momentum()-decay[0]->momentum();
  vmom.rescaleMass();
  for(hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(ix=decay.size();ix>1;--ix) ihel[ix]=(hhel%constants[ix-1])/constants[ix];
    vWave=VectorWaveFunction(vmom,vec,hadron[hhel]*UnitRemoval::InvE,outgoing);
    for(unsigned int if1 = 0; if1 < 2; ++if1) {
      for(unsigned int if2 = 0; if2 < 2; ++if2) {
	ihel[0]=if1;
	ihel[1]=if2;
	if(!ferm) swap(ihel[0],ihel[1]);
	(*newME)(ihel) = FFVPtr_->evaluate(scale,wave_[if1],wavebar_[if2],vWave);
      }
    }
  }
  // store the matrix element
  ME(newME);
  // multiply by the CKM element
  int iq,ia;
  weakCurrent()->decayModeInfo(mode(),iq,ia);
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);
  }
  pre /= 0.125*sqr(FFVPtr_->weakCoupling(scale));
  double output(0.5*pre*ckm*(ME()->contract(rho_)).real()*
		sqr(SM().fermiConstant()*UnitRemoval::E2));
  return output;
}
 
Energy FFVCurrentDecayer::partialWidth(tPDPtr inpart, tPDPtr outa,
				       vector<tPDPtr> currout) {
  vector<long> id;
  id.push_back(inpart->id());
  id.push_back(outa->id());
  for(unsigned int ix=0;ix<currout.size();++ix) id.push_back(currout[ix]->id());
  bool cc;
  int mode=modeNumber(cc,id);
  imode(mode);
  return initializePhaseSpaceMode(mode,true,true);  
}
