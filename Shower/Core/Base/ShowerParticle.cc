// -*- C++ -*-
//
// ShowerParticle.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerParticle class.
//

#include "ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include <ostream>

using namespace Herwig;

PPtr ShowerParticle::clone() const {
  return new_ptr(*this);
}

PPtr ShowerParticle::fullclone() const {
  return new_ptr(*this);
}

ClassDescription<ShowerParticle> ShowerParticle::initShowerParticle;
// Definition of the static class description member.

void ShowerParticle::vetoEmission(ShowerPartnerType, Energy scale) {
  scales_.QED         = min(scale,scales_.QED        );
  scales_.QED_noAO    = min(scale,scales_.QED_noAO   );
  scales_.QCD_c       = min(scale,scales_.QCD_c      );
  scales_.QCD_c_noAO  = min(scale,scales_.QCD_c_noAO );
  scales_.QCD_ac      = min(scale,scales_.QCD_ac     );
  scales_.QCD_ac_noAO = min(scale,scales_.QCD_ac_noAO);
}

void ShowerParticle::addPartner(EvolutionPartner in ) {
  partners_.push_back(in); 
}

namespace {

LorentzRotation boostToShower(Lorentz5Momentum & porig, tShowerBasisPtr basis) {
  LorentzRotation output; 
  assert(basis);
  if(basis->frame()==ShowerBasis::BackToBack) {
    // we are doing the evolution in the back-to-back frame
    // work out the boostvector
    Boost boostv(-(basis->pVector()+basis->nVector()).boostVector());
    // momentum of the parton
    Lorentz5Momentum ptest(basis->pVector());
    // construct the Lorentz boost
    output = LorentzRotation(boostv);
    ptest *= output;
    Axis axis(ptest.vect().unit());
    // now rotate so along the z axis as needed for the splitting functions
    if(axis.perp2()>1e-10) {
      double sinth(sqrt(1.-sqr(axis.z())));
      output.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    }
    else if(axis.z()<0.) {
      output.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    porig = output*basis->pVector();
    porig.setX(ZERO);
    porig.setY(ZERO);
  }
  else {
    output = LorentzRotation(-basis->pVector().boostVector());
    porig = output*basis->pVector();
    porig.setX(ZERO);
    porig.setY(ZERO);
    porig.setZ(ZERO);
  }
  return output;
}

RhoDMatrix bosonMapping(ShowerParticle & particle,
			const Lorentz5Momentum & porig,
			VectorSpinPtr vspin,
			const LorentzRotation & rot) {
  // rotate the original basis
  vector<LorentzPolarizationVector> sbasis;
  for(unsigned int ix=0;ix<3;++ix) {
    sbasis.push_back(vspin->getProductionBasisState(ix));
    sbasis.back().transform(rot);
  }
  // splitting basis
  vector<LorentzPolarizationVector> fbasis;
  bool massless(particle.id()==ParticleID::g||particle.id()==ParticleID::gamma);
  VectorWaveFunction wave(porig,particle.dataPtr(),outgoing);
  for(unsigned int ix=0;ix<3;++ix) {
    if(massless&&ix==1) {
      fbasis.push_back(LorentzPolarizationVector());
    }
    else {
      wave.reset(ix);
      fbasis.push_back(wave.wave());
    }
  }
  // work out the mapping
  RhoDMatrix mapping=RhoDMatrix(PDT::Spin1,false);
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      mapping(ix,iy)= sbasis[iy].dot(fbasis[ix].conjugate());
      if(particle.id()<0)
	mapping(ix,iy)=conj(mapping(ix,iy));
    }
  }
  // \todo need to fix this
  mapping = RhoDMatrix(PDT::Spin1,false);
  if(massless) {
    mapping(0,0) = 1.;
    mapping(2,2) = 1.;
  }
  else {
    mapping(0,0) = 1.;
    mapping(1,1) = 1.;
    mapping(2,2) = 1.;
  }
  return mapping;
}

RhoDMatrix fermionMapping(ShowerParticle & particle,
			  const Lorentz5Momentum & porig,
			  FermionSpinPtr fspin,
			  const LorentzRotation & rot) {
  // extract the original basis states
  vector<LorentzSpinor<SqrtEnergy> > sbasis;
  for(unsigned int ix=0;ix<2;++ix) {
    sbasis.push_back(fspin->getProductionBasisState(ix));
    sbasis.back().transform(rot);
  }
  // calculate the states in the splitting basis
  vector<LorentzSpinor<SqrtEnergy> > fbasis;
  SpinorWaveFunction wave(porig,particle.dataPtr(),
			  particle.id()>0 ? incoming : outgoing);
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    fbasis.push_back(wave.dimensionedWave());
  }
  RhoDMatrix mapping=RhoDMatrix(PDT::Spin1Half,false);
  for(unsigned int ix=0;ix<2;++ix) {
    if(fbasis[0].s2()==complex<SqrtEnergy>()) {
      mapping(ix,0) = sbasis[ix].s3()/fbasis[0].s3();
      mapping(ix,1) = sbasis[ix].s2()/fbasis[1].s2();
    }
    else {
      mapping(ix,0) = sbasis[ix].s2()/fbasis[0].s2();
      mapping(ix,1) = sbasis[ix].s3()/fbasis[1].s3();
    }
  }
  return mapping;
}

FermionSpinPtr createFermionSpinInfo(ShowerParticle & particle,
				     const Lorentz5Momentum & porig,
				     const LorentzRotation & rot,
				     Helicity::Direction dir) {
  // calculate the splitting basis for the branching
  // and rotate back to construct the basis states
  LorentzRotation rinv = rot.inverse();
  SpinorWaveFunction wave;
  if(particle.id()>0)
    wave=SpinorWaveFunction(porig,particle.dataPtr(),incoming);
  else
    wave=SpinorWaveFunction(porig,particle.dataPtr(),outgoing);
  FermionSpinPtr fspin = new_ptr(FermionSpinInfo(particle.momentum(),dir==outgoing));
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    LorentzSpinor<SqrtEnergy> basis = wave.dimensionedWave();
    basis.transform(rinv);
    fspin->setBasisState(ix,basis);
    fspin->setDecayState(ix,basis);
  }
  particle.spinInfo(fspin);
  return fspin;
}

VectorSpinPtr createVectorSpinInfo(ShowerParticle & particle,
				   const Lorentz5Momentum & porig,
				   const LorentzRotation & rot,
				   Helicity::Direction dir) {
  // calculate the splitting basis for the branching
  // and rotate back to construct the basis states
  LorentzRotation rinv = rot.inverse();
  bool massless(particle.id()==ParticleID::g||particle.id()==ParticleID::gamma);
  VectorWaveFunction wave(porig,particle.dataPtr(),dir);
  VectorSpinPtr vspin = new_ptr(VectorSpinInfo(particle.momentum(),dir==outgoing));
  for(unsigned int ix=0;ix<3;++ix) {
    LorentzPolarizationVector basis;
    if(massless&&ix==1) {
      basis = LorentzPolarizationVector();
    }
    else {
      wave.reset(ix);
      basis = wave.wave();
    }
    basis *= rinv;
    vspin->setBasisState(ix,basis);
    vspin->setDecayState(ix,basis);
  }
  particle.spinInfo(vspin);
  vspin->  DMatrix() = RhoDMatrix(PDT::Spin1);
  vspin->rhoMatrix() = RhoDMatrix(PDT::Spin1);
  if(massless) {
    vspin->  DMatrix()(0,0) = 0.5;
    vspin->rhoMatrix()(0,0) = 0.5;
    vspin->  DMatrix()(2,2) = 0.5;
    vspin->rhoMatrix()(2,2) = 0.5;
  }
  return vspin;
}

}

RhoDMatrix ShowerParticle::extractRhoMatrix(bool forward) {
  // get the spin density matrix and the mapping
  RhoDMatrix mapping;
  SpinPtr inspin;
  bool needMapping = getMapping(inspin,mapping);
  // set the decayed flag
  inspin->decay();
  // get the spin density matrix
  RhoDMatrix rho = forward ? inspin->rhoMatrix() : inspin->DMatrix();
  // map to the shower basis if needed  
  if(needMapping) {
    RhoDMatrix rhop(rho.iSpin(),false);
    for(int ixb=0;ixb<rho.iSpin();++ixb) {
  	for(int iyb=0;iyb<rho.iSpin();++iyb) {
  	  if(mapping(iyb,ixb)==0.)continue;
	  for(int iya=0;iya<rho.iSpin();++iya) {
            if(rho(iya,iyb)==0.)continue;
	    for(int ixa=0;ixa<rho.iSpin();++ixa) {
  	      rhop(ixa,ixb) += rho(iya,iyb)*mapping(iya,ixa)*conj(mapping(iyb,ixb));
  	    }
  	  }
  	}
    }
    rhop.normalize();
    rho = rhop;
  }
  return rho;
}

bool ShowerParticle::getMapping(SpinPtr & output, RhoDMatrix & mapping) {
  // if the particle is not from the hard process
  if(!this->perturbative()) {
    // mapping is the identity
    output=this->spinInfo();
    mapping=RhoDMatrix(this->dataPtr()->iSpin());
    if(output) {
      return false;
    }
    else {
      Lorentz5Momentum porig;
      LorentzRotation rot = boostToShower(porig,showerBasis());
      Helicity::Direction dir = this->isFinalState() ? outgoing : incoming;
      if(this->dataPtr()->iSpin()==PDT::Spin0) {
	assert(false);
      }
      else if(this->dataPtr()->iSpin()==PDT::Spin1Half) {
	output = createFermionSpinInfo(*this,porig,rot,dir);
      }
      else if(this->dataPtr()->iSpin()==PDT::Spin1) {
	output = createVectorSpinInfo(*this,porig,rot,dir);
      }
      else {
	assert(false);
      }
      return false;
    }
  }
  // if particle is final-state and is from the hard process
  else if(this->isFinalState()) {
    assert(this->perturbative()==1 || this->perturbative()==2);
    // get transform to shower frame
    Lorentz5Momentum porig;
    LorentzRotation rot = boostToShower(porig,showerBasis());
    // the rest depends on the spin of the particle
    PDT::Spin spin(this->dataPtr()->iSpin());
    mapping=RhoDMatrix(spin,false);
    // do the spin dependent bit
    if(spin==PDT::Spin0) {
      ScalarSpinPtr sspin=dynamic_ptr_cast<ScalarSpinPtr>(this->spinInfo());
      if(!sspin) {
	ScalarWaveFunction::constructSpinInfo(this,outgoing,true);
      }
      output=this->spinInfo();
      return false;
    }
    else if(spin==PDT::Spin1Half) {
      FermionSpinPtr fspin=dynamic_ptr_cast<FermionSpinPtr>(this->spinInfo());
      // spin info exists get information from it
      if(fspin) {
	output=fspin;
	mapping = fermionMapping(*this,porig,fspin,rot);
	return true;
      }
      // spin info does not exist create it
      else {
	output = createFermionSpinInfo(*this,porig,rot,outgoing);
	return false;
      }
    }
    else if(spin==PDT::Spin1) {
      VectorSpinPtr vspin=dynamic_ptr_cast<VectorSpinPtr>(this->spinInfo());
      // spin info exists get information from it
      if(vspin) {
	output=vspin;
	mapping = bosonMapping(*this,porig,vspin,rot);
	return true; 
      }
      else {
	output = createVectorSpinInfo(*this,porig,rot,outgoing);
	return false;
      }
    }
    // not scalar/fermion/vector
    else
      assert(false);
  }
  // incoming to hard process
  else if(this->perturbative()==1 && !this->isFinalState()) {
    // get the basis vectors
    // get transform to shower frame
    Lorentz5Momentum porig;
    LorentzRotation rot = boostToShower(porig,showerBasis());
    porig *= this->x();
    // the rest depends on the spin of the particle
    PDT::Spin spin(this->dataPtr()->iSpin());
    mapping=RhoDMatrix(spin);
    // do the spin dependent bit
    if(spin==PDT::Spin0) {
      cerr << "testing spin 0 not yet implemented " << endl;
      assert(false);
    }
    // spin-1/2
    else if(spin==PDT::Spin1Half) {
      FermionSpinPtr fspin=dynamic_ptr_cast<FermionSpinPtr>(this->spinInfo());
      // spin info exists get information from it
      if(fspin) {
	output=fspin;
	mapping = fermionMapping(*this,porig,fspin,rot);
	return true;
      }
      // spin info does not exist create it
      else {
	output = createFermionSpinInfo(*this,porig,rot,incoming);
	return false;
      }
    }
    // spin-1
    else if(spin==PDT::Spin1) {
      VectorSpinPtr vspin=dynamic_ptr_cast<VectorSpinPtr>(this->spinInfo());
      // spinInfo exists map it
      if(vspin) {
	output=vspin;
	mapping = bosonMapping(*this,porig,vspin,rot);
	return true;
      }
      // create the spininfo
      else {
	output = createVectorSpinInfo(*this,porig,rot,incoming);
	return false;
      }
    }
    assert(false);
  }
  // incoming to decay
  else if(this->perturbative() == 2 && !this->isFinalState()) { 
    // get the basis vectors
    Lorentz5Momentum porig;
    LorentzRotation rot=boostToShower(porig,showerBasis());
    // the rest depends on the spin of the particle
    PDT::Spin spin(this->dataPtr()->iSpin());
    mapping=RhoDMatrix(spin);
    // do the spin dependent bit
    if(spin==PDT::Spin0) {
      cerr << "testing spin 0 not yet implemented " << endl;
      assert(false);
    }
    // spin-1/2
    else if(spin==PDT::Spin1Half) {
    //   FermionSpinPtr fspin=dynamic_ptr_cast<FermionSpinPtr>(this->spinInfo());
    //   // spin info exists get information from it
    //   if(fspin) {
    // 	output=fspin;
    // 	mapping = fermionMapping(*this,porig,fspin,rot);
    // 	return true;
    //   // spin info does not exist create it
    //   else {
    // 	output = createFermionSpinInfo(*this,porig,rot,incoming);
    // 	return false;
    //   }
    // }
      assert(false);
    }
    // // spin-1
    // else if(spin==PDT::Spin1) {
    //   VectorSpinPtr vspin=dynamic_ptr_cast<VectorSpinPtr>(this->spinInfo());
    //   // spinInfo exists map it
    //   if(vspin) {
    // 	output=vspin;
    // 	mapping = bosonMapping(*this,porig,vspin,rot);
    // 	return true;
    //   }
    //   // create the spininfo
    //   else {
    // 	output = createVectorSpinInfo(*this,porig,rot,incoming);
    // 	return false;
    //   }
    // }
    // assert(false);
    assert(false);
  }
  else
    assert(false);
  return true;
}

void ShowerParticle::constructSpinInfo(bool timeLike) {
  // now construct the required spininfo and calculate the basis states
  PDT::Spin spin(dataPtr()->iSpin());
  if(spin==PDT::Spin0) {
    ScalarWaveFunction::constructSpinInfo(this,outgoing,timeLike);
  }
  // calculate the basis states and construct the SpinInfo for a spin-1/2 particle
  else if(spin==PDT::Spin1Half) {
    // outgoing particle
    if(id()>0) {
      vector<LorentzSpinorBar<SqrtEnergy> > stemp;
      SpinorBarWaveFunction::calculateWaveFunctions(stemp,this,outgoing);
      SpinorBarWaveFunction::constructSpinInfo(stemp,this,outgoing,timeLike);
    }
    // outgoing antiparticle
    else {
      vector<LorentzSpinor<SqrtEnergy> > stemp;
      SpinorWaveFunction::calculateWaveFunctions(stemp,this,outgoing);
      SpinorWaveFunction::constructSpinInfo(stemp,this,outgoing,timeLike);
    }
  }
  // calculate the basis states and construct the SpinInfo for a spin-1 particle
  else if(spin==PDT::Spin1) {
    bool massless(id()==ParticleID::g||id()==ParticleID::gamma);
    vector<Helicity::LorentzPolarizationVector> vtemp;
    VectorWaveFunction::calculateWaveFunctions(vtemp,this,outgoing,massless);
    VectorWaveFunction::constructSpinInfo(vtemp,this,outgoing,timeLike,massless);
  }
  else {
    throw Exception() << "Spins higher than 1 are not yet implemented in " 
		      << "FS_QtildaShowerKinematics1to2::constructVertex() "
		      << Exception::runerror;
  }
}

void ShowerParticle::initializeDecay() {
  Lorentz5Momentum p, n, ppartner, pcm;
  assert(perturbative()!=1);
  // this is for the initial decaying particle
  if(perturbative()==2) {
    ShowerBasisPtr newBasis;
    p = momentum();
    Lorentz5Momentum ppartner(partner()->momentum());
    // removed to make inverse recon work properly
    //if(partner->thePEGBase()) ppartner=partner->thePEGBase()->momentum();
    pcm=ppartner;
    Boost boost(p.findBoostToCM());
    pcm.boost(boost);
    n = Lorentz5Momentum( ZERO,0.5*p.mass()*pcm.vect().unit()); 
    n.boost( -boost);
    newBasis = new_ptr(ShowerBasis());
    newBasis->setBasis(p,n,ShowerBasis::Rest);
    showerBasis(newBasis,false);
  }
  else {
    showerBasis(dynamic_ptr_cast<ShowerParticlePtr>(parents()[0])->showerBasis(),true);
  }
}

void ShowerParticle::initializeInitialState(PPtr parent) {
  // For the time being we are considering only 1->2 branching
  Lorentz5Momentum p, n, pthis, pcm;
  assert(perturbative()!=2);
  if(perturbative()==1) {
    ShowerBasisPtr newBasis;
    // find the partner and its momentum
    if(!partner()) return;
    if(partner()->isFinalState()) {
      Lorentz5Momentum pa = -momentum()+partner()->momentum();
      Lorentz5Momentum pb =  momentum();
      Energy scale=parent->momentum().t();
      Lorentz5Momentum pbasis(ZERO,parent->momentum().vect().unit()*scale);
      Axis axis(pa.vect().unit());
      LorentzRotation rot;
      double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
      if(axis.perp2()>1e-20) {
	rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
	rot.rotateX(Constants::pi);
      }
      if(abs(1.-pa.e()/pa.vect().mag())>1e-6) rot.boostZ( pa.e()/pa.vect().mag());
      pb *= rot;
      if(pb.perp2()/GeV2>1e-20) {
	Boost trans = -1./pb.e()*pb.vect();
	trans.setZ(0.);
	rot.boost(trans);
      }
      pbasis *=rot;
      rot.invert();
      n = rot*Lorentz5Momentum(ZERO,-pbasis.vect());
      p = rot*Lorentz5Momentum(ZERO, pbasis.vect());
    }
    else {
      pcm = parent->momentum();
      p = Lorentz5Momentum(ZERO, pcm.vect());
      n = Lorentz5Momentum(ZERO, -pcm.vect());
    }
    newBasis = new_ptr(ShowerBasis());
    newBasis->setBasis(p,n,ShowerBasis::BackToBack);
    showerBasis(newBasis,false);
  } 
  else {
    showerBasis(dynamic_ptr_cast<ShowerParticlePtr>(children()[0])->showerBasis(),true);
  }
}

void ShowerParticle::initializeFinalState() {
  // set the basis vectors
  Lorentz5Momentum p,n;
  if(perturbative()!=0) {
    ShowerBasisPtr newBasis;
    // find the partner() and its momentum
    if(!partner()) return;
    Lorentz5Momentum ppartner(partner()->momentum());
    // momentum of the emitting particle
    p = momentum();
    Lorentz5Momentum pcm;
    // if the partner is a final-state particle then the reference
    // vector is along the partner in the rest frame of the pair
    if(partner()->isFinalState()) {
      Boost boost=(p + ppartner).findBoostToCM();
      pcm = ppartner;
      pcm.boost(boost);
      n = Lorentz5Momentum(ZERO,pcm.vect());
      n.boost( -boost);
    }
    else if(!partner()->isFinalState()) {
      // if the partner is an initial-state particle then the reference
      // vector is along the partner which should be massless
      if(perturbative()==1) {
	n = Lorentz5Momentum(ZERO,ppartner.vect());
      }
      // if the partner is an initial-state decaying particle then the reference
      // vector is along the backwards direction in rest frame of decaying particle
      else {
	Boost boost=ppartner.findBoostToCM();
	pcm = p;
	pcm.boost(boost);
	n = Lorentz5Momentum( ZERO, -pcm.vect()); 
	n.boost( -boost);
      } 
    }
    newBasis = new_ptr(ShowerBasis());
    newBasis->setBasis(p,n,ShowerBasis::BackToBack);
    showerBasis(newBasis,false);
  }
  else if(initiatesTLS()) {
    showerBasis(dynamic_ptr_cast<ShowerParticlePtr>(parents()[0]->children()[0])->showerBasis(),true);
  }
  else  {
    showerBasis(dynamic_ptr_cast<ShowerParticlePtr>(parents()[0])->showerBasis(),true);
  }
}

void ShowerParticle::setShowerMomentum(bool timeLike) {
  Energy m = this->mass() > ZERO ? this->mass() : this->data().mass();
  // calculate the momentum of the assuming on-shell
  Energy2 pt2 = sqr(this->showerParameters().pt);
  double alpha = timeLike ? this->showerParameters().alpha : this->x();
  double beta = 0.5*(sqr(m) + pt2 - sqr(alpha)*showerBasis()->pVector().m2())/(alpha*showerBasis()->p_dot_n());
  Lorentz5Momentum porig=showerBasis()->sudakov2Momentum(alpha,beta,
							 this->showerParameters().ptx,
							 this->showerParameters().pty);
  porig.setMass(m);
  this->set5Momentum(porig);
}
