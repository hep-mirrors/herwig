// -*- C++ -*-
//
// SMWDecayer.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWDecayer class.
//

#include "SMWDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMWDecayer::SMWDecayer() 
  : quarkWeight_(6,0.), leptonWeight_(3,0.) {
  quarkWeight_[0]  = 1.01596;
  quarkWeight_[1]  = 0.0537308;
  quarkWeight_[2]  = 0.0538085;
  quarkWeight_[3]  = 1.01377;
  quarkWeight_[4]  = 1.45763e-05;
  quarkWeight_[5]  = 0.0018143;
  leptonWeight_[0] = 0.356594;
  leptonWeight_[1] = 0.356593;
  leptonWeight_[2] = 0.356333;
  // intermediates
  generateIntermediates(false);
}

void SMWDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig++ StandardModel object in"
				  << "SMWDecayer::doinit()"
				  << Exception::runerror;
  FFWvertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFW());
  // make sure they are initialized
  FFWvertex_->init();
  // now set up the decay modes
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  vector<double> wgt(0);
  // W modes
  extpart[0]=getParticleData(ParticleID::Wplus);
  // loop for the quarks
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(!FFWvertex_->allowed(-ix,iy,ParticleID::Wminus))
	throw InitException() << "SMWDecayer::doinit() the W vertex" 
			      << "cannot handle all the quark modes" 
			      << Exception::abortnow;
      extpart[1] = getParticleData(-ix);
      extpart[2] = getParticleData( iy);
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,quarkWeight_[iz],wgt);
      ++iz;
    }
  }
  // loop for the leptons
  for(int ix=11;ix<17;ix+=2) {
    // check that the combination of particles is allowed
    if(!FFWvertex_->allowed(-ix,ix+1,ParticleID::Wminus))
      throw InitException() << "SMWDecayer::doinit() the W vertex" 
			    << "cannot handle all the lepton modes" 
			    << Exception::abortnow;
    extpart[1] = getParticleData(-ix);
    extpart[2] = getParticleData(ix+1);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    addMode(mode,leptonWeight_[(ix-11)/2],wgt);
  }
}

int SMWDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(abs(id0)!=ParticleID::Wplus) return imode;
  int idd(0),idu(0);
  if(abs(id1)%2==1&&abs(id2)%2==0) {
    idd=abs(id1);
    idu=abs(id2);
  }
  else if(abs(id1)%2==0&&abs(id2)%2==1) {
    idd=abs(id2);
    idu=abs(id1);
  }
  if(idd==0&&idu==0) {
    return imode;
  }
  else if(idd<=5) {
    imode=idd+idu/2-2;
  }
  else {
    imode=(idd-1)/2+1;
  }
  cc= (id0==ParticleID::Wminus);
  return imode;
}

void SMWDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFWvertex_ 
     << quarkWeight_ << leptonWeight_;
  
}

void SMWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWvertex_ 
     >> quarkWeight_ >> leptonWeight_;
}

ClassDescription<SMWDecayer> SMWDecayer::initSMWDecayer;
// Definition of the static class description member.

void SMWDecayer::Init() {

  static ClassDocumentation<SMWDecayer> documentation
    ("The SMWDecayer class is the implementation of the decay"
     " of the W boson to the Standard Model fermions.");
  static ParVector<SMWDecayer,double> interfaceWquarkMax
    ("QuarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWDecayer::quarkWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWDecayer,double> interfaceWleptonMax
    ("LeptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWDecayer::leptonWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

}


// return the matrix element squared
double SMWDecayer::me2(const int, const Particle & inpart,
			const ParticleVector & decay,
			MEOption meopt) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    ME(DecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half));
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
    return 0.;
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(_wavebar,decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(_wave   ,decay[ianti],outgoing);
  // compute the matrix element
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	if(iferm>ianti) ME()(vhel,ia,ifm)=
	  FFWvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
	else            ME()(vhel,ifm,ia)=
	  FFWvertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME().contract(_rho)).real()*UnitRemoval::E2/scale;
  if(abs(decay[0]->id())<=6) output*=3.;
  if(decay[0]->hasColour())      decay[0]->antiColourNeighbour(decay[1]);
  else if(decay[1]->hasColour()) decay[1]->antiColourNeighbour(decay[0]);
  return output;
}

void SMWDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<6) quarkWeight_ [ix]=mode(ix)->maxWeight();
      else     leptonWeight_[ix-6]=mode(ix)->maxWeight();
    }
  }
}

void SMWDecayer::dataBaseOutput(ofstream & output,
				 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<quarkWeight_.size();++ix) {
    output << "set " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "set " << name() << ":LeptonMax " << ix << " "
	   << leptonWeight_[ix] << "\n";
  }
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
