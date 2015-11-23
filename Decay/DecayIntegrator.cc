// -*- C++ -*-
//
// DecayIntegrator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayIntegrator class.
//
// Author: Peter Richardson
// 

#include "DecayIntegrator.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Utilities/Debug.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "DecayPhaseSpaceMode.h"
#include "Herwig/PDT/WidthCalculatorBase.h"
#include "ThePEG/Interface/Reference.h"

using namespace Herwig; 

ParticleVector DecayIntegrator::decay(const Particle & parent,
				      const tPDVector & children) const {
  // return empty vector if products heavier than parent
  Energy mout(ZERO);
  for(tPDVector::const_iterator it=children.begin();
      it!=children.end();++it) mout+=(**it).massMin();
  if(mout>parent.mass()) return ParticleVector();
  // generate the decay
  bool cc;
  _imode = modeNumber(cc,parent.dataPtr(),children);
  return _modes[_imode]->generate(_generateinter,cc,parent);
}
  
void DecayIntegrator::persistentOutput(PersistentOStream & os) const {
  os << _modes << _niter << _npoint << _ntry << _photongen << _generateinter;
}
  
void DecayIntegrator::persistentInput(PersistentIStream & is, int) {
  is >> _modes >> _niter >> _npoint >> _ntry >> _photongen >> _generateinter;
}
  
AbstractClassDescription<DecayIntegrator> DecayIntegrator::initDecayIntegrator;
// Definition of the static class description member.
  
void DecayIntegrator::Init() {
    
  static RefVector<DecayIntegrator,DecayPhaseSpaceMode> interfaceModes
    ("Modes",
     "The phase space integration modes.",
     &DecayIntegrator::_modes, 0, false, false, true, true); 
  
  static ClassDocumentation<DecayIntegrator> documentation
    ("The DecayIntegrator class is a base decayer class "
     "including a multi-channel integrator.");
  
  static Parameter<DecayIntegrator,int> interfaceIteration
    ("Iteration",
     "Number of iterations for the initialization of the phase space",
     &DecayIntegrator::_niter, 10, 0, 100,
     false, false, true);  
  
  static Parameter<DecayIntegrator,int> interfacePoints
    ("Points",
     "number of phase space points to generate in the initialisation.",
     &DecayIntegrator::_npoint, 10000, 1, 1000000000,
     false, false, true);
  
  static Parameter<DecayIntegrator,int> interfaceNtry
    ("Ntry",
     "Number of attempts to generate the decay",
     &DecayIntegrator::_ntry, 500, 0, 100000,
     false, false, true);

  static Reference<DecayIntegrator,DecayRadiationGenerator> interfacePhotonGenerator
    ("PhotonGenerator",
     "Object responsible for generating photons in the decay.",
     &DecayIntegrator::_photongen, false, false, true, true, false);
 
  static Switch<DecayIntegrator,bool> interfaceGenerateIntermediates
    ("GenerateIntermediates",
     "Whether or not to include intermediate particles in the output",
     &DecayIntegrator::_generateinter, false, false, false);
  static SwitchOption interfaceGenerateIntermediatesNoIntermediates
    (interfaceGenerateIntermediates,
     "No",
     "Don't include the intermediates",
     false);
  static SwitchOption interfaceGenerateIntermediatesIncludeIntermediates
    (interfaceGenerateIntermediates,
     "Yes",
     "include the intermediates",
     true);
}

// output info on the integrator
ostream & Herwig::operator<<(ostream & os, const DecayIntegrator & decay) {
  os << "The integrator has " << decay._modes.size() << " modes"  << endl;
  for(unsigned int ix=0;ix<decay._modes.size();++ix) {
    os << "Information on mode " << ix << endl;
    os << *(decay._modes[ix]);
  }
  return os;
}

// generate the momenta for the decay
ParticleVector DecayIntegrator::generate(bool inter,bool cc,
					 const unsigned int & imode,
					 const Particle & inpart) const {
  _imode=imode;
  return _modes[imode]->generate(inter,cc,inpart);
}  

// initialization for a run
void DecayIntegrator::doinitrun() {
  HwDecayerBase::doinitrun();
  if ( initialize() && Debug::level > 1 ) 
    CurrentGenerator::current().log() << "Start of the initialisation for " 
				      << this->name() << "\n";
  for(unsigned int ix=0;ix<_modes.size();++ix) {
    if(!_modes[ix]) continue;
    _modes[ix]->initrun();
    _imode=ix;
    _modes[ix]->initializePhaseSpace(initialize());
  }
  // CurrentGenerator::current().log() << *this << "\n";
}

// add a new mode
void DecayIntegrator::addMode(DecayPhaseSpaceModePtr in,double maxwgt,
				     const vector<double> inwgt) const {
  _modes.push_back(in);
  if(in) {
    in->setMaxWeight(maxwgt);
    in->setWeights(inwgt);
    in->setIntegrate(_niter,_npoint,_ntry);
    in->init();
  }
}

// reset the properities of all intermediates
void DecayIntegrator::resetIntermediate(tcPDPtr part, Energy mass, Energy width) {
  if(!part) return;
  for(unsigned int ix=0,N=_modes.size();ix<N;++ix) {
    _modes[ix]->resetIntermediate(part,mass,width);
  }
} 

WidthCalculatorBasePtr 
DecayIntegrator::threeBodyMEIntegrator(const DecayMode &) const {
  return WidthCalculatorBasePtr();
}

// set the code for the partial width
void DecayIntegrator::setPartialWidth(const DecayMode & dm, int imode) {
  vector<int> extid;
  tcPDPtr cc,cc2;
  int nfound(0),ifound,nmax(1),id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  if(dm.parent()->CC()) nmax=2;
  if(_modes.size()==0) return;
  do {
    if(!_modes[ix]) {
      ++ix;
      continue;
    }
    cc = _modes[ix]->externalParticles(0)->CC();
    tmax = cc ? 1 : 2;
    for(iz=0;iz<tmax;++iz) {
      ifound=-1;
      extid.clear();
      // check the parent
      if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==0) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  extid.push_back(_modes[ix]->externalParticles(iy)->id());
	}
      }
      else if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==1) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  cc2=_modes[ix]->externalParticles(iy)->CC();
	  extid.push_back( cc2 ? cc2->id() : _modes[ix]->externalParticles(iy)->id());
	}
      }
      else if(cc&&dm.parent()->id()==cc->id()) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  cc2 = _modes[ix]->externalParticles(iy)->CC();
	  extid.push_back( cc2 ? cc2->id() : _modes[ix]->externalParticles(iy)->id());
	}
      }
      // if the parents match
      if(!extid.empty()) {
	vector<bool> matched(extid.size(),false);
	bool done;
	nmatched=0;
	ParticleMSet::const_iterator pit = dm.products().begin();
	do {
	  id=(**pit).id();
	  done=false;
	  iy=1;
	  do {
	    if(id==extid[iy]&&!matched[iy]) {
	      matched[iy]=true;
	      ++nmatched;
	      done=true;
	    }
	    ++iy;
	  }
	  while(iy<extid.size()&&!done);
	  ++pit;
	}
	while(pit!=dm.products().end());
	if(nmatched==extid.size()-1) {
	  ifound=ix;
	  ++nfound;
	}
      }
      if(ifound>=0) _modes[ifound]->setPartialWidth(imode);
    }
    ++ix;
  }
  while(nfound<nmax&&ix<_modes.size());
}

int DecayIntegrator::findMode(const DecayMode & dm) {
  int imode(-1);
  vector<int> extid;
  tcPDPtr cc,cc2;
  bool found(false);
  int id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  if(_modes.size()==0) return -1;
  do {
    if(!_modes[ix]) {
      ++ix;
      continue;
    }
    cc = _modes[ix]->externalParticles(0)->CC();
    tmax=1;if(!cc){++tmax;}
    for(iz=0;iz<tmax;++iz) {
      extid.clear();
      // check the parent
      if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==0) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  extid.push_back(_modes[ix]->externalParticles(iy)->id());
	}
      }
      else if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==1) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  cc2=_modes[ix]->externalParticles(iy)->CC();
	  extid.push_back( cc2 ? cc2->id() : _modes[ix]->externalParticles(iy)->id());
	}
      }
      else if(cc&&dm.parent()->id()==cc->id()) {
	for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy) {
	  cc = _modes[ix]->externalParticles(iy)->CC();
	  extid.push_back( cc ? cc->id() : _modes[ix]->externalParticles(iy)->id());
	}
      }
      // if the parents match
      if(!extid.empty()) {
	vector<bool> matched(extid.size(),false);
	bool done;
	nmatched=0;
	ParticleMSet::const_iterator pit = dm.products().begin();
	do {
	  id=(**pit).id();
	  done=false;
	  iy=1;
	  do {
	    if(id==extid[iy]&&!matched[iy]) {
	      matched[iy]=true;
	      ++nmatched;
	      done=true;
	    }
	    ++iy;
	  }
	  while(iy<extid.size()&&!done);
	  ++pit;
	}
	while(pit!=dm.products().end());
	if(nmatched==extid.size()-1) {
	  imode=ix;
	  found=true;
	}
      }
    }
    ++ix;
  }
  while(!found&&ix<_modes.size());
  return imode;
}

// output the information for the database
void DecayIntegrator::dataBaseOutput(ofstream & output,bool header) const {
  // header for MySQL
  if(header) output << "update decayers set parameters=\"";
  output << "newdef " << name() << ":Iteration " << _niter << "\n";
  output << "newdef " << name() << ":Ntry " << _ntry << "\n";
  output << "newdef " << name() << ":Points " << _npoint << "\n";
  //if(_photongen){;}
  output << "newdef " << name() << ":GenerateIntermediates " << _generateinter << " \n";
  // footer for MySQL
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
  }
}
  
// pointer to a mode
tDecayPhaseSpaceModePtr DecayIntegrator::mode(unsigned int ix) {
  return _modes[ix];
}

tcDecayPhaseSpaceModePtr DecayIntegrator::mode(unsigned int ix) const {
  return _modes[ix];
}
 
Energy DecayIntegrator::initializePhaseSpaceMode(unsigned int imode,bool init, bool onShell) const{
  tcDecayPhaseSpaceModePtr cmodeptr=mode(imode);
  tDecayPhaseSpaceModePtr modeptr = const_ptr_cast<tDecayPhaseSpaceModePtr>(cmodeptr);
  modeptr->init();
  return modeptr->initializePhaseSpace(init,onShell);
}

// the matrix element to be integrated for the me
double DecayIntegrator::threeBodyMatrixElement(const int,const Energy2,
					       const Energy2,
					       const Energy2,const Energy2,
					       const Energy, const Energy, 
					       const Energy) const {
  throw DecayIntegratorError() 
    << "Calling the virtual DecayIntegrator::threeBodyMatrixElement"
    << "method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

// the differential three body decay rate with one integral performed
InvEnergy DecayIntegrator::threeBodydGammads(const int, const Energy2,
					     const Energy2,
					     const Energy, const Energy, 
					     const Energy) const {
  throw DecayIntegratorError() 
    << "Calling the virtual DecayIntegrator::threeBodydGammads()" 
    <<"method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

double DecayIntegrator::oneLoopVirtualME(unsigned int ,
					 const Particle &, 
					 const ParticleVector &) {
  throw DecayIntegratorError()
    << "DecayIntegrator::oneLoopVirtualME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}

InvEnergy2 DecayIntegrator::realEmissionME(unsigned int,
					   const Particle &, 
					   ParticleVector &,
					   unsigned int,
					   double, double, 
					   const LorentzRotation &,
					   const LorentzRotation &) {
  throw DecayIntegratorError()
    << "DecayIntegrator::realEmmisionME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}
