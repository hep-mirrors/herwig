// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VertexBase class.
//

#include "VertexBase.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

namespace Herwig{
namespace Helicity{
using namespace ThePEG;
using namespace ThePEG::Helicity;
        
void VertexBase::doinit() throw(InitException) {
  Interfaced::doinit();
  // get the particle data points for the external particles
  tPDPtr pin;
  if(_npoint>=3)
    {
      for(unsigned int ix=0;ix<size();++ix)
	{
	  pin=getParticleData(_iparticlea[ix]);
	  if(!pin) throw InitException() 
	    << "Unknown particle ID=" << _iparticlea[ix] 
	    << "requested in VertexBase::doinit() for " << this->fullName() 
	    << Exception::runerror;
	  _particlea.push_back(pin);
	  pin=getParticleData(_iparticleb[ix]);
	  if(!pin) throw InitException() 
	    << "Unknown particle ID=" << _iparticleb[ix] 
	    << "requested in VertexBase::doinit() for " << this->fullName() 
	    << Exception::runerror;
	  _particleb.push_back(pin);
	  pin=getParticleData(_iparticlec[ix]);
	  if(!pin) throw InitException() 
	    << "Unknown particle ID=" << _iparticlec[ix] 
	    << "requested in VertexBase::doinit() for " << this->fullName() 
	    << Exception::runerror;
	  _particlec.push_back(pin);
	}
    }
  if(_npoint>=4)
    {
      for(unsigned int ix=0;ix<size();++ix)
	{
	  pin=getParticleData(_iparticled[ix]);
	  if(!pin) throw InitException() 
	    << "Unknown particle ID=" << _iparticled[ix] 
	    << "requested in VertexBase::doinit() for " << this->fullName() 
	    << Exception::runerror;
	  _particled.push_back(pin);
	}
    }
  if(_npoint==5)
    {
      for(unsigned int ix=0;ix<size();++ix)
	{
	  pin=getParticleData(_iparticlee[ix]);
	  if(!pin) throw InitException() 
	    << "Unknown particle ID=" << _iparticlee[ix] 
	    << "requested in VertexBase::doinit() for " << this->fullName() 
	    << Exception::runerror;
	  _particlee.push_back(pin);
	}
    }
  // set up the incoming and outgoing particles
  setIncoming();
  setOutgoing();
}
    
void VertexBase::persistentOutput(PersistentOStream & os) const {
  os << _npoint << _nsize << _ispin << _inpart << _iinpart << _outpart << _ioutpart 
     << _iparticlea << _iparticleb << _iparticlec << _iparticled << _iparticlee
     << _particlea << _particleb << _particlec << _particled << _particlee 
     << _calckinematics;
}

void VertexBase::persistentInput(PersistentIStream & is, int) {
  is >> _npoint >> _nsize >> _ispin >> _inpart >> _iinpart >> _outpart >> _ioutpart 
     >> _iparticlea >> _iparticleb >> _iparticlec >> _iparticled >> _iparticlee
     >> _particlea >> _particleb >> _particlec >> _particled >> _particlee
     >> _calckinematics;
}
    
AbstractClassDescription<VertexBase> VertexBase::initVertexBase;
// Definition of the static class description member.
  
void VertexBase::Init() {
 
  static Parameter<VertexBase,unsigned int> interfacenpoint
    ("NPoint",
     "The number of extermal particles interacting at the Vertex.",
     &VertexBase::_npoint, 3, 3, 5, false, false, true);
  
  static Switch<VertexBase,bool> interfaceCalculateKinematics
    ("CalculateKinematics",
     "Calculate kinematic invariants at the vertices. This is"
     " mainly needed for loop vertices.",
     &VertexBase::_calckinematics, false, false, false);
  static SwitchOption interfaceCalculateKinematicsCalculate
    (interfaceCalculateKinematics,
     "Calculate",
     "Calculate the kinematics",
     true);
  static SwitchOption interfaceCalculateKinematicsNoKinematics
    (interfaceCalculateKinematics,
     "NoKinematics",
     "Do not calculate the kinematics",
     false);

  static ParVector<VertexBase,int> interfaceispin
    ("Spin",
     "The spins of the external particles interacting at the Vertex.",
     &VertexBase::_ispin,
     0, 0, 0, 0, 5, false, false, true);
  
  static ParVector<VertexBase,int> interfacefirstparticle
    ("FirstParticle",
     "Possible first particles for the Vertex",
     &VertexBase::_iparticlea,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VertexBase,int> interfacesecondparticle
    ("SecondParticle",
     "Possible second particles for the Vertex",
     &VertexBase::_iparticleb,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VertexBase,int> interfacethirdparticle
    ("ThirdParticle",
     "Possible third particles for the Vertex",
     &VertexBase::_iparticlec,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VertexBase,int> interfacefourthparticle
    ("FourthParticle",
     "Possible fourth particles for the Vertex",
     &VertexBase::_iparticled,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VertexBase,int> interfacefifthparticle
    ("FifthParticle",
     "Possible fifth particles for the Vertex",
     &VertexBase::_iparticlee,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ClassDocumentation<VertexBase> documentation
    ("The VertexBase class is designed to be the base class"
     "of all vertices in Herwig++");

}

// find particles with a given id    
vector<PDPtr> VertexBase::search(unsigned int iloc,int idd)
{
  vector<PDPtr> out;
  if(iloc>=_npoint)
    {throw HelicityConsistencyError() << "VertexBase::search Invalid _particle "
				      << "index for ilist search" 
				      << Exception::abortnow;}
  else
    {
      bool found;
      for(unsigned int ix=0;ix<_particlea.size();++ix)
	{
	  found=false;
	  switch (iloc)
	    {
	    case 0:
	      if(_particlea[ix]->id()==idd){found=true;}
	      break;
	    case 1:
	      if(_particleb[ix]->id()==idd){found=true;}
	      break;
	    case 2:
	      if(_particlec[ix]->id()==idd){found=true;}
	      break;
	    case 3:
	      if(_particled[ix]->id()==idd){found=true;}
	      break;
	    case 4:
	      if(_particlee[ix]->id()==idd){found=true;}
	      break;
	    default:
	      std::cerr << "Invalid _particle index for ilist search" << std::endl;
	      break;
	    }
	  if(found)
	    {
	      if(_npoint>=3)
		{
		  out.push_back(_particlea[ix]);
		  out.push_back(_particleb[ix]);
		  out.push_back(_particlec[ix]);
		}
	      if(_npoint>=4){out.push_back(_particled[ix]);}
	      if(_npoint==5){out.push_back(_particlee[ix]);}
	    }
	}
    }
  return out;
}

// check a given combination is allowed for a three point vertex
bool VertexBase::allowed(int ida, int idb, int idc)
{
  if(_npoint!=3)
    {
      throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				       << " a three point Vertex" << Exception::warning;
      return false;
    }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0){return false;}
  else
    {
      int iloop=out.size()/_npoint, iy;
      for(int ix=0;ix<iloop;++ix)
	{
	  iy = ix*_npoint;
	  if((out[iy+1])->id()==idb && (out[iy+2])->id()==idc){return true;}
	}
    }
  return false;
}

// check a given combination is allowed for a four point vertex
bool VertexBase::allowed(int ida, int idb, int idc, int idd)
{
  if(_npoint!=4)
    {
      throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				       << " a four point Vertex" << Exception::warning;
      return false;
    }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0)
    {return false;}
  else
    {
      int iloop=out.size()/_npoint, iy;
      for(int ix=0;ix<iloop;++ix)
	{
	  iy = ix*_npoint;
	  if(out[iy+1]->id()==idb && out[iy+2]->id()==idc && out[iy+3]->id()==idd)
	    {return true;}
	}
    }
  return false;
}

bool VertexBase::allowed(int ida, int idb, int idc, int idd, int ide)
{
  if(_npoint!=5)
    {
      throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				       << " a five point Vertex" << Exception::warning;
      return false;
    }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0)
    {return false;}
  else
    {
      int iloop=out.size()/_npoint, iy;
      for(int ix=0;ix<iloop;++ix)
	{
	  iy = ix*_npoint;
	  if(out[iy+1]->id()==idb && out[iy+2]->id()==idc && 
	     out[iy+3]->id()==idd && out[iy+4]->id()==ide)
	    {return true;}
	}
    }
  return false;
}
  
// output the information
ostream & operator<<(ostream & os, const VertexBase & in)
{
  os << "Information on Vertex" << endl;
  os << "This is an " << in._npoint << " vertex" << endl;
  if(in._calckinematics){os << "The kinematic invariants are calculated" << endl;}
  else{os << "The kinematics invariants are not calculated" << endl;}
  os << " Particles allowed for this Vertex" << endl;
  for(unsigned int ix=0;ix<in._particlea.size();++ix)
    {
      if(in._npoint==3)
	{
	  os << in._particlea[ix]->id() << "   "
	     << in._particleb[ix]->id() << "   "
	     << in._particlec[ix]->id() << "   " << endl;
	}
      else if(in._npoint==4)
	{
	  os << in._particlea[ix]->id() << "   "
	     << in._particleb[ix]->id() << "   "
	     << in._particlec[ix]->id() << "   "
	     << in._particled[ix]->id() << "   " << endl;
	}
      else if(in._npoint==5)
	{
	  os << in._particlea[ix]->id() << "   "
	     << in._particleb[ix]->id() << "   "
	     << in._particlec[ix]->id() << "   "
	     << in._particled[ix]->id() << "   "
	     << in._particlee[ix]->id() << "   " << endl;
	}
      else
	{
	  os << "Invalid number of external particles" << endl;
	}
    }
  return os;
}


// add particle to the list for a three point vertex
void VertexBase::add(int ia ,int ib ,int ic)
{
  if(_npoint==3)
    {
      // add to the PDG code lists
      _iparticlea.push_back(ia);_iparticleb.push_back(ib);_iparticlec.push_back(ic);
      // add to the Particle data pointer lists
      tPDPtr pin;
      pin=getParticleData(ia);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ia 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlea.push_back(pin);
      pin=getParticleData(ib);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ib 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particleb.push_back(pin);
      pin=getParticleData(ic);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ic 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlec.push_back(pin);
    }
  else
    {throw HelicityConsistencyError() << "This is a " << _npoint 
				      << " vertex cannot add three particles" 
				      << Exception::abortnow;}
  // add to the list of outgoing particles
  if(!outgoing(ia)){_outpart.push_back(_particlea[_nsize]);_ioutpart.push_back(ia);}
  if(!outgoing(ib)){_outpart.push_back(_particleb[_nsize]);_ioutpart.push_back(ib);}
  if(!outgoing(ic)){_outpart.push_back(_particlec[_nsize]);_ioutpart.push_back(ic);}
  // add to the list of incoming particles  
  if(_particlea[_nsize]->CC())
    {
      if(!incoming(-ia))
	{
	  _inpart.push_back(_particlea[_nsize]->CC());
	  _iinpart.push_back(-ia);
	}
    }
  else
    {
      if(!incoming(ia))
	{
	  _inpart.push_back(_particlea[_nsize]);
	  _iinpart.push_back(ia);
	}
    }
  if(_particleb[_nsize]->CC())
    {
      if(!incoming(-ib))
	{
	  _inpart.push_back(_particleb[_nsize]->CC());
	  _iinpart.push_back(-ib);
	}
    }
  else
    {
      if(!incoming(ib))
	{
	  _inpart.push_back(_particleb[_nsize]);
	  _iinpart.push_back(ib);
	}
    }
  if(_particlec[_nsize]->CC())
    {
      if(!incoming(-ic))
	{
	  _inpart.push_back(_particlec[_nsize]->CC());
	  _iinpart.push_back(-ic);
	}
    }
  else
    {
      if(!incoming(ic))
	{
	  _inpart.push_back(_particlec[_nsize]);
	  _iinpart.push_back(ic);
	}
    }
  // increment the size of the arrays
  ++_nsize;
}

// add particle to the list for a four point vertex
void VertexBase::add(int ia,int ib,int ic,int id)
{
  if(_npoint==4)
    {
      // add to the PDG code lists
      _iparticlea.push_back(ia);_iparticleb.push_back(ib);
      _iparticlec.push_back(ic);_iparticled.push_back(id);
      // add to the Particle data pointer lists
      tPDPtr pin;
      pin=getParticleData(ia);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ia 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlea.push_back(pin);
      pin=getParticleData(ib);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ib 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particleb.push_back(pin);
      pin=getParticleData(ic);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ic 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlec.push_back(pin);
      pin=getParticleData(id);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << id 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particled.push_back(pin);
    }
  {throw HelicityConsistencyError() << "This is a " << _npoint 
				    << " vertex cannot add four particles" 
				    << Exception::abortnow;}
  // add to the list of outgoing particles
  if(!outgoing(ia)){_outpart.push_back(_particlea[_nsize]);_ioutpart.push_back(ia);}
  if(!outgoing(ib)){_outpart.push_back(_particleb[_nsize]);_ioutpart.push_back(ib);}
  if(!outgoing(ic)){_outpart.push_back(_particlec[_nsize]);_ioutpart.push_back(ic);}
  if(!outgoing(id)){_outpart.push_back(_particled[_nsize]);_ioutpart.push_back(id);}
  // add to the list of incoming particles  
  if(_particlea[_nsize]->CC())
    {
      if(!incoming(-ia))
	{
	  _inpart.push_back(_particlea[_nsize]->CC());
	  _iinpart.push_back(-ia);
	}
    }
  else
    {
      if(!incoming(ia))
	{
	  _inpart.push_back(_particlea[_nsize]);
	  _iinpart.push_back(ia);
	}
    }
  if(_particleb[_nsize]->CC())
    {
      if(!incoming(-ib))
	{
	  _inpart.push_back(_particleb[_nsize]->CC());
	  _iinpart.push_back(-ib);
	}
    }
  else
    {
      if(!incoming(ib))
	{
	  _inpart.push_back(_particleb[_nsize]);
	  _iinpart.push_back(ib);
	}
    }
  if(_particlec[_nsize]->CC())
    {
      if(!incoming(-ic))
	{
	  _inpart.push_back(_particlec[_nsize]->CC());
	  _iinpart.push_back(-ic);
	}
    }
  else
    {
      if(!incoming(ic))
	{
	  _inpart.push_back(_particlec[_nsize]);
	  _iinpart.push_back(ic);
	}
    }
  if(_particled[_nsize]->CC())
    {
      if(!incoming(-id))
	{
	  _inpart.push_back(_particled[_nsize]->CC());
	_iinpart.push_back(-id);
	}
    }
  else
    {
      if(!incoming(id))
	{
	  _inpart.push_back(_particled[_nsize]);
	_iinpart.push_back(id);
	}	
    }
  // increment the size of the arrays
  ++_nsize;
}

// add particle to the list for a five point vertex
void VertexBase::add(int ia,int ib,int ic,int id, int ie)
{
  if(_npoint==5)
    {
      // add to the PDG code lists
      _iparticlea.push_back(ia);_iparticleb.push_back(ib);
      _iparticlec.push_back(ic);_iparticled.push_back(id);_iparticlee.push_back(ie);
      // add to the Particle data pointer lists
      tPDPtr pin;
      pin=getParticleData(ia);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ia 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlea.push_back(pin);
      pin=getParticleData(ib);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ib 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particleb.push_back(pin);
      pin=getParticleData(ic);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ic 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlec.push_back(pin);
      pin=getParticleData(id);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << id 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particled.push_back(pin);
      pin=getParticleData(ie);
      if(!pin) throw InitException() 
	<< "Unknown particle ID=" << ie 
	<< "requested in VertexBase::add() for " << this->fullName() 
	<< Exception::runerror;
      _particlee.push_back(pin);
    }
  {throw HelicityConsistencyError() << "This is a " << _npoint 
				    << " vertex cannot add five particles" 
				    << Exception::abortnow;}
  // add to the list of outgoing particles
  if(!outgoing(ia)){_outpart.push_back(_particlea[_nsize]);_ioutpart.push_back(ia);}
  if(!outgoing(ib)){_outpart.push_back(_particleb[_nsize]);_ioutpart.push_back(ib);}
  if(!outgoing(ic)){_outpart.push_back(_particlec[_nsize]);_ioutpart.push_back(ic);}
  if(!outgoing(id)){_outpart.push_back(_particled[_nsize]);_ioutpart.push_back(id);}
  if(!outgoing(ie)){_outpart.push_back(_particlee[_nsize]);_ioutpart.push_back(ie);}
  // add to the list of incoming particles  
  if(_particlea[_nsize]->CC())
    {
      if(!incoming(-ia))
	{
	  _inpart.push_back(_particlea[_nsize]->CC());
	  _iinpart.push_back(-ia);
	}
    }
  else
    {
      if(!incoming(ia))
	{
	  _inpart.push_back(_particlea[_nsize]);
	  _iinpart.push_back(ia);
	}
    }
  if(_particleb[_nsize]->CC())
    {
      if(!incoming(-ib))
	{
	  _inpart.push_back(_particleb[_nsize]->CC());
	  _iinpart.push_back(-ib);
	}
    }
  else
    {
      if(!incoming(ib))
	{
	  _inpart.push_back(_particleb[_nsize]);
	  _iinpart.push_back(ib);
	}
    }
  if(_particlec[_nsize]->CC())
    {
      if(!incoming(-ic))
	{
	  _inpart.push_back(_particlec[_nsize]->CC());
	  _iinpart.push_back(-ic);
	}
    }
  else
    {
      if(!incoming(ic))
	{
	  _inpart.push_back(_particlec[_nsize]);
	  _iinpart.push_back(ic);
	}
    }
  if(_particled[_nsize]->CC())
    {
      if(!incoming(-id))
	{
	  _inpart.push_back(_particled[_nsize]->CC());
	_iinpart.push_back(-id);
	}
    }
  else
    {
      if(!incoming(id))
	{
	  _inpart.push_back(_particled[_nsize]);
	_iinpart.push_back(id);
	}	
    }
  if(_particlee[_nsize]->CC())
    {
      if(!incoming(-ie))
	{
	  _inpart.push_back(_particlee[_nsize]->CC());
	  _iinpart.push_back(-ie);
	}
    }
  else
    {
      if(!incoming(ie))
	{
	  _inpart.push_back(_particlee[_nsize]);
	_iinpart.push_back(ie);
	}
    }
  // increment the size of the arrays
  ++_nsize;
}


// set the list of outgoing particles
void VertexBase::setOutgoing()
{
  if(_outpart.size()==0)
    {
      if(_npoint>=3)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      if(!outgoing(_iparticlea[iy]))
		{
		  _outpart.push_back(_particlea[iy]);
		  _ioutpart.push_back(_iparticlea[iy]);
		}
	      if(!outgoing(_iparticleb[iy]))
		{
		  _outpart.push_back(_particleb[iy]);
		  _ioutpart.push_back(_iparticleb[iy]);
		}
	      if(!outgoing(_iparticlec[iy]))
		{
		  _outpart.push_back(_particlec[iy]);
		  _ioutpart.push_back(_iparticlec[iy]);
		}
	    }
	}
      if(_npoint>=4)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      if(!outgoing(_iparticled[iy]))
		{
		  _outpart.push_back(_particled[iy]);
		  _ioutpart.push_back(_iparticled[iy]);
		}
	    }
	}
      if(_npoint==5)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      if(!outgoing(_iparticlee[iy]))
		{
		  _outpart.push_back(_particlee[iy]);
		  _ioutpart.push_back(_iparticlee[iy]);
		}
	    }
	}
    }
  else
    {throw HelicityConsistencyError() << "VertexBase::setOutgoing " 
				      << "Outgoing particles already set" 
				      << Exception::abortnow;} 
}

// set the list of incoming particles
void VertexBase::setIncoming()
{
  if(_inpart.size()==0)
    {
      PDPtr temp;
      if(_npoint>=3)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      temp =_particlea[iy]->CC();
	      if(temp)
		{
		  if(!incoming(-_iparticlea[iy]))
		    {
		      _inpart.push_back(temp);
		      _iinpart.push_back(-_iparticlea[iy]);
		    }
		}
	      else
		{
		  if(!incoming(_iparticlea[iy]))
		    {
		      _inpart.push_back(_particlea[iy]);
		      _iinpart.push_back(_iparticlea[iy]);
		    }
		}
	      temp =_particleb[iy]->CC();
	      if(temp)
		{
		  if(!incoming(-_iparticleb[iy]))
		    {
		      _inpart.push_back(temp);
		      _iinpart.push_back(-_iparticleb[iy]);
		    }
		}
	      else
		{
		  if(!incoming(_iparticleb[iy]))
		    {
		      _inpart.push_back(_particleb[iy]);
		      _iinpart.push_back(_iparticleb[iy]);
		    }
		}
	      temp =_particlec[iy]->CC();
	      if(temp)
		{
		  if(!incoming(-_iparticlec[iy]))
		    {
		      _inpart.push_back(temp);
		      _iinpart.push_back(-_iparticlec[iy]);
		    }
		}
	      else
		{
		  if(!incoming(_iparticlec[iy]))
		    {
		      _inpart.push_back(_particlec[iy]);
		      _iinpart.push_back(_iparticlec[iy]);
		    }
		}
	    }

	}
      if(_npoint>=4)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      temp =_particled[iy]->CC();
	      if(temp)
		{
		  if(!incoming(-_iparticled[iy]))
		    {
		      _inpart.push_back(temp);
		      _iinpart.push_back(-_iparticled[iy]);
		    }
		}
	      else
		{
		  if(!incoming(_iparticled[iy]))
		    {
		      _inpart.push_back(_particled[iy]);
		      _iinpart.push_back(_iparticled[iy]);
		    }
		}
	    }
	}
      if(_npoint==5)
	{
	  for(unsigned int iy=0;iy<size();++iy)
	    {
	      temp =_particlee[iy]->CC();
	      if(temp)
		{
		  if(!incoming(-_iparticlee[iy]))
		    {
		      _inpart.push_back(temp);
		      _iinpart.push_back(-_iparticlee[iy]);
		    }
		}
	      else
		{
		  if(!incoming(_iparticlee[iy]))
		    {
		      _inpart.push_back(_particlee[iy]);
		      _iinpart.push_back(_iparticlee[iy]);
		    }
		}
	    }
	}
    }   
  else{throw HelicityConsistencyError() << "VertexBase::setIncoming " 
					<< "Outgoing particles already set" 
					<< Exception::abortnow;}   
}
 
}
}

