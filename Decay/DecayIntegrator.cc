// -*- C++ -*-
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
#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/PDT/GenericMassGenerator.h"
#include "ThePEG/Utilities/Timer.h"
#include "DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/WidthCalculatorBase.h"

namespace Herwig {
using namespace ThePEG;

// copy constructor
DecayIntegrator::DecayIntegrator(const DecayIntegrator & x)
  : Decayer(x),_niter(x._niter),_npoint(x._npoint), _ntry(x._ntry),
    _Initialize(x._Initialize), _modes(x._modes) {}

// destructor
DecayIntegrator::~DecayIntegrator() {}

// dummy accept method
bool DecayIntegrator::accept(const DecayMode & dm) const {return false;}

// dummy decay method
ParticleVector DecayIntegrator::decay(const DecayMode & dm,
				      const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}
  
void DecayIntegrator::persistentOutput(PersistentOStream & os) const {
  os <<_Initialize << _modes << _niter << _npoint << _ntry;
  }
  
void DecayIntegrator::persistentInput(PersistentIStream & is, int) {
  is >>_Initialize >> _modes >> _niter >> _npoint >> _ntry;
}
  
AbstractClassDescription<DecayIntegrator> DecayIntegrator::initDecayIntegrator;
// Definition of the static class description member.
  
void DecayIntegrator::Init() {
    
  static RefVector<DecayIntegrator,DecayPhaseSpaceMode> interfaceModes
    ("Modes",
     "The phase space integration modes.",
     &DecayIntegrator::_modes, 0, false, false, true, true); 
  
  static ClassDocumentation<DecayIntegrator> documentation
    ("The \\classname{DecayIntegrator} class is a base decayer class "
     "including a multi-channel integrator.");
  
  static Switch<DecayIntegrator,bool> interfaceInitialize
    ("Initialize",
     "Initialization of the phase space calculation",
     &DecayIntegrator::_Initialize, false, false, false);

  static SwitchOption interfaceInitializeon
    (interfaceInitialize,
     "on",
     "At initialisation find max weight and optimise the integration",
     true);

  static SwitchOption interfaceInitializeoff
    (interfaceInitialize,
     "off",
     "Use the maximum weight and channel weights supplied for the integration",
     false);

  static Parameter<DecayIntegrator,int> interfaceIteration
    ("Iteration",
     "Number of iterations for the initialization of the phase space",
     &DecayIntegrator::_niter, 10, 0, 100,
     false, false, true);  
  
  static Parameter<DecayIntegrator,int> interfacePoints
    ("Points",
     "number of phase space points to generate in the initialisation.",
     &DecayIntegrator::_npoint, 10000, 100, 100000000,
     false, false, true);
  
  static Parameter<DecayIntegrator,int> interfaceNtry
    ("Ntry",
     "Number of attempts to generate the decay",
     &DecayIntegrator::_ntry, 500, 0, 100000,
     false, false, true);
   
}

// output info on the integrator
ostream & operator<<(ostream & os, const DecayIntegrator & decay)
{
  os << "The integrator has " << decay._modes.size() << " modes"  << endl;
 for(unsigned int ix=0;ix<decay._modes.size();++ix)
   {
     os << "Information on mode " << ix << endl;
     os << *(decay._modes[ix]);
   }
 return os;
}

// generate the momenta for the decay
ParticleVector DecayIntegrator::generate(bool inter,bool cc, const unsigned int & imode,
					 const Particle & inpart) const
 {
   _imode=imode;
   return _modes[imode]->generate(inter,cc,inpart);
 }  


  // initialization for a run
void DecayIntegrator::doinitrun() {
  Decayer::doinitrun();
  CurrentGenerator::current().log() << "testing start of the initialisation for " 
				    << this->fullName() << endl;
  for(unsigned int ix=0;ix<_modes.size();++ix)
    {
      _modes[ix]->initrun();
      _imode=ix;_modes[ix]->initializePhaseSpace(_Initialize);
    }
}

// add a new mode
void DecayIntegrator::addMode(DecayPhaseSpaceModePtr in,double maxwgt,
				     const vector<double> inwgt) const
{
  _modes.push_back(in);
  in->setMaxWeight(maxwgt);
  in->setWeights(inwgt);
  in->setIntegrate(_niter,_npoint,_ntry);
  in->init();
}

// reset the properities of all intermediates
void DecayIntegrator::resetIntermediate(tcPDPtr part, Energy mass, Energy width)
{
  for(unsigned int ix=0,N=_modes.size();ix<N;++ix)
    {_modes[ix]->resetIntermediate(part,mass,width);}
} 

bool DecayIntegrator::twoBodyMEcode(const DecayMode & dm, int &, double &) const
{
  throw DecayIntegratorError() << "Calling the virtual DecayIntegrator::twoBodyMEcode"
			       << " method this must be overwritten in the classes "
			       << "inheriting from DecayIntegrator where it is needed"
			       << Exception::runerror;
}

// the matrix element to be integrated for the me
double DecayIntegrator::threeBodyMatrixElement(int,Energy2,Energy2,Energy2,Energy2,
					       Energy,Energy,Energy)
{
  throw DecayIntegratorError() << "Calling the virtual DecayIntegrator::threeBodyMatrixElement"
			       << "method. This must be overwritten in the classes "
			       << "inheriting from DecayIntegrator where it is needed"
			       << Exception::runerror;
}

  // the differential three body decay rate with one integral performed
double DecayIntegrator::threeBodydGammads(int,Energy2,Energy2,Energy,Energy,Energy)
{
  throw DecayIntegratorError() << "Calling the virtual DecayIntegrator::threeBodydGammads()" 
			       <<"method. This must be overwritten in the classes "
			       << "inheriting from DecayIntegrator where it is needed"
			       << Exception::runerror;
}

WidthCalculatorBasePtr DecayIntegrator::threeBodyMEIntegrator(const DecayMode & dm) const
{return WidthCalculatorBasePtr();}


// set the code for the partial width
void DecayIntegrator::setPartialWidth(const DecayMode & dm, int imode)
{
  vector<int> extid;
  tcPDPtr cc,cc2;
  int nfound(0),ifound,nmax(1),id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  if(dm.parent()->CC()){nmax=2;}
  if(_modes.size()==0){return;}
  do 
    {
      cc = _modes[ix]->externalParticles(0)->CC();
      tmax=1;if(!cc){++tmax;}
      for(iz=0;iz<tmax;++iz)
	{
	  ifound=-1;
	  extid.resize(0);
	  // check the parent
	  if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==0)
	    {for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{extid.push_back(_modes[ix]->externalParticles(iy)->id());}}
	  else if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==1)
	    {
	      for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{
		  cc2=_modes[ix]->externalParticles(iy)->CC();
		  if(cc2){extid.push_back(cc2->id());}
		  else{extid.push_back(_modes[ix]->externalParticles(iy)->id());}
		}
	    }
	  else if(cc&&dm.parent()->id()==cc->id())
	    {
	      for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{
		  cc = _modes[ix]->externalParticles(iy)->CC();
		  if(cc){extid.push_back(cc->id());}
		  else{extid.push_back(_modes[ix]->externalParticles(iy)->id());}
		}
	    }
	  // if the parents match
	  if(!extid.empty())
	    {
	      vector<bool> matched(extid.size(),false);bool done;
	      nmatched=0;
	      ParticleMSet::const_iterator pit = dm.products().begin();
	      do
		{
		  id=(**pit).id();
		  done=false;
		  iy=1;
		  do 
		    {
		      if(id==extid[iy]&&!matched[iy])
			{matched[iy]=true;++nmatched;done=true;}
		      ++iy;
		    }
		  while(iy<extid.size()&&!done);
		  ++pit;
		}
	      while(pit!=dm.products().end());
	      if(nmatched==extid.size()-1){ifound=ix;++nfound;}
	    }
	  if(ifound>=0){_modes[ifound]->setPartialWidth(imode);}
	}
      ++ix;
    }
  while(nfound<nmax&&ix<_modes.size());
}

int DecayIntegrator::findMode(const DecayMode & dm)
{
  int imode(-1);
  vector<int> extid;
  tcPDPtr cc,cc2;
  bool found(false);
  int id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  cout << "testing the number of modes " << _modes.size() << endl;

  if(_modes.size()==0){return -1;}
  do 
    {
      cc = _modes[ix]->externalParticles(0)->CC();
      tmax=1;if(!cc){++tmax;}
      for(iz=0;iz<tmax;++iz)
	{
	  extid.resize(0);
	  // check the parent
	  if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==0)
	    {for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{extid.push_back(_modes[ix]->externalParticles(iy)->id());}}
	  else if(dm.parent()->id()==_modes[ix]->externalParticles(0)->id()&&iz==1)
	    {
	      for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{
		  cc2=_modes[ix]->externalParticles(iy)->CC();
		  if(cc2){extid.push_back(cc2->id());}
		  else{extid.push_back(_modes[ix]->externalParticles(iy)->id());}
		}
	    }
	  else if(cc&&dm.parent()->id()==cc->id())
	    {
	      for(iy=0,N=_modes[ix]->numberofParticles();iy<N;++iy)
		{
		  cc = _modes[ix]->externalParticles(iy)->CC();
		  if(cc){extid.push_back(cc->id());}
		  else{extid.push_back(_modes[ix]->externalParticles(iy)->id());}
		}
	    }
	  // if the parents match
	  if(!extid.empty())
	    {
	      vector<bool> matched(extid.size(),false);bool done;
	      nmatched=0;
	      ParticleMSet::const_iterator pit = dm.products().begin();
	      do
		{
		  id=(**pit).id();
		  done=false;
		  iy=1;
		  do 
		    {
		      if(id==extid[iy]&&!matched[iy])
			{matched[iy]=true;++nmatched;done=true;}
		      ++iy;
		    }
		  while(iy<extid.size()&&!done);
		  ++pit;
		}
	      while(pit!=dm.products().end());
	      if(nmatched==extid.size()-1){imode=ix;found=true;}
	    }
	}
      ++ix;
    }
  while(!found&&ix<_modes.size());
  cout << "testing found mode " << _modes[imode]->externalParticles(0)->PDGName() 
	<< " -> ";
  for(ix=1;ix<_modes[imode]->numberofParticles();++ix)
    {cout << _modes[imode]->externalParticles(ix)->PDGName() << " ";}
  cout << endl;
  return imode;
}

// output the information for the database
void DecayIntegrator::dataBaseOutput(ofstream & output)
{output << " where ThePEGName=\" " << fullName() << "\";";}

// pointer to a mode
tDecayPhaseSpaceModePtr DecayIntegrator::mode(unsigned int ix){return _modes[ix];}
tcDecayPhaseSpaceModePtr DecayIntegrator::mode(unsigned int ix) const
{return _modes[ix];}
}
