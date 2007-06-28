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

using namespace Herwig;
using namespace Herwig::Helicity;
using namespace ThePEG;
using namespace ThePEG::Helicity;
        
void VertexBase::doinit() throw(InitException) {
  Interfaced::doinit();
  _nsize=0;
  // get the particle data points for the external particles
  tPDPtr pin;
  if(_npoint==3) {
    for(unsigned int ix=0;ix<_iparticlea.size();++ix) {
      tPDPtr pin[3] = {getParticleData(_iparticlea[ix]),
		       getParticleData(_iparticleb[ix]),
		       getParticleData(_iparticlec[ix])};
      if(!pin[0]||!pin[1]||!pin[2]) continue;
      _particles.push_back(vector<PDPtr>(pin,pin+3));
      ++_nsize;
    }
  }
  if(_npoint>=4) {
    for(unsigned int ix=0;ix<_iparticlea.size();++ix) {
      tPDPtr pin[4] = {getParticleData(_iparticlea[ix]),
		       getParticleData(_iparticleb[ix]),
		       getParticleData(_iparticlec[ix]),
		       getParticleData(_iparticled[ix])};
      if(!pin[0]||!pin[1]||!pin[2]||!pin[3]) continue;
      _particles.push_back(vector<PDPtr>(pin,pin+4));
      ++_nsize;
    }
  }
  if(_npoint==5) {
    for(unsigned int ix=0;ix<_iparticlea.size();++ix) {
      tPDPtr pin[5] = {getParticleData(_iparticlea[ix]),
		       getParticleData(_iparticleb[ix]),
		       getParticleData(_iparticlec[ix]),
		       getParticleData(_iparticled[ix]),
		       getParticleData(_iparticlee[ix])};
      if(!pin[0]||!pin[1]||!pin[2]||!pin[3]||!pin[4]) continue;
      _particles.push_back(vector<PDPtr>(pin,pin+5));
      ++_nsize;
    }
  }
  // set up the incoming and outgoing particles
  setIncoming();
  setOutgoing();
}
    
void VertexBase::persistentOutput(PersistentOStream & os) const {
  os << _npoint << _nsize << _ispin << _inpart << _iinpart << _outpart << _ioutpart 
     << _iparticlea << _iparticleb << _iparticlec << _iparticled << _iparticlee
     << _particles << _calckinematics;
}

void VertexBase::persistentInput(PersistentIStream & is, int) {
  is >> _npoint >> _nsize >> _ispin >> _inpart >> _iinpart >> _outpart >> _ioutpart 
     >> _iparticlea >> _iparticleb >> _iparticlec >> _iparticled >> _iparticlee
     >> _particles >> _calckinematics;
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
vector<PDPtr> VertexBase::search(unsigned int iloc,int idd) {
  vector<PDPtr> out;
  if(iloc>=_npoint) 
    throw HelicityConsistencyError() << "VertexBase::search Invalid _particle "
				     << "index for ilist search" 
				     << Exception::abortnow;
  bool found;
  for(unsigned int ix=0;ix<_particles.size();++ix) {
    found=_particles[ix][iloc]->id()==idd;
    if(found) {
      for(unsigned int iy=0;iy<_particles[ix].size();++iy) {
	out.push_back(_particles[ix][iy]);
      }
    }
  }
  return out;
}

// check a given combination is allowed for a three point vertex
bool VertexBase::allowed(int ida, int idb, int idc) {
  if(_npoint!=3) {
    throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				     << " a three point Vertex" << Exception::warning;
    return false;
  }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0) return false;
  int iloop=out.size()/_npoint, iy;
  for(int ix=0;ix<iloop;++ix) {
    iy = ix*_npoint;
    if((out[iy+1])->id()==idb && (out[iy+2])->id()==idc) return true;
  }
  return false;
}

// check a given combination is allowed for a four point vertex
bool VertexBase::allowed(int ida, int idb, int idc, int idd) {
  if(_npoint!=4) {
    throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				     << " a four point Vertex" << Exception::warning;
    return false;
  }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0) return false;
  int iloop=out.size()/_npoint, iy;
  for(int ix=0;ix<iloop;++ix) {
    iy = ix*_npoint;
    if(out[iy+1]->id()==idb && out[iy+2]->id()==idc && out[iy+3]->id()==idd)
      return true;
  }
  return false;
}

bool VertexBase::allowed(int ida, int idb, int idc, int idd, int ide) {
  if(_npoint!=5) {
    throw HelicityConsistencyError() << "VertexBase::allowed Not allowed as not"
				     << " a five point Vertex" << Exception::warning;
    return false;
  }
  vector<PDPtr> out = search(0,ida);
  if(out.size()==0) return false;
  int iloop=out.size()/_npoint, iy;
  for(int ix=0;ix<iloop;++ix) {
    iy = ix*_npoint;
    if(out[iy+1]->id()==idb && out[iy+2]->id()==idc && 
       out[iy+3]->id()==idd && out[iy+4]->id()==ide)
      return true;
  }
  return false;
}
  
// output the information
ostream & Herwig::Helicity::operator<<(ostream & os, const VertexBase & in) {
  os << "Information on Vertex" << endl;
  os << "This is an " << in._npoint << " vertex" << endl;
  if(in._calckinematics){os << "The kinematic invariants are calculated" << endl;}
  else{os << "The kinematics invariants are not calculated" << endl;}
  os << " Particles allowed for this Vertex" << endl;
  for(unsigned int ix=0;ix<in._particles.size();++ix) {
    for(unsigned int iy=0;iy<in._particles[ix].size();++iy) {
      os << in._particles[ix][iy]->id() << "   ";
    }
    os << "\n";
  }
  return os;
}


// add particle to the list for a three point vertex
void VertexBase::add(int ia ,int ib ,int ic) {
  if(_npoint!=3) throw HelicityConsistencyError() << "This is a " << _npoint 
				      << " vertex cannot add three particles" 
				      << Exception::abortnow;
  // add to the PDG code lists
  _iparticlea.push_back(ia);
  _iparticleb.push_back(ib);
  _iparticlec.push_back(ic);
  tPDPtr pin[3]={getParticleData(ia),
		 getParticleData(ib),
		 getParticleData(ic)};
  if(!pin[0]||!pin[1]||!pin[2]) return;
  // add to the list of outgoing particles
  _particles.push_back(vector<PDPtr>(pin,pin+3));
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(!outgoing(_particles.back()[ix]->id())) {
      _outpart.push_back(_particles.back()[ix]);
      _ioutpart.push_back(_particles.back()[ix]->id());
    }
  }
  // add to the list of incoming particles  
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(_particles.back()[ix]->CC()) {
      if(!incoming(_particles.back()[ix]->CC()->id())) {
	_inpart.push_back(_particles.back()[ix]->CC());
	_iinpart.push_back(_particles.back()[ix]->CC()->id());
      }
    }
    else {
      if(!incoming(_particles.back()[ix]->id())) {
	_inpart.push_back(_particles.back()[ix]);
	_iinpart.push_back(_particles.back()[ix]->id());
      }
    }
  }
  // increment the size of the arrays
  ++_nsize;
}

// add particle to the list for a four point vertex
void VertexBase::add(int ia,int ib,int ic,int id) {
  if(_npoint!=4)
    throw HelicityConsistencyError() << "This is a " << _npoint 
				     << " vertex cannot add four particles" 
				     << Exception::abortnow;
  // add to the PDG code lists
  _iparticlea.push_back(ia);
  _iparticleb.push_back(ib);
  _iparticlec.push_back(ic);
  _iparticled.push_back(id);
  tPDPtr pin[4]={getParticleData(ia),
		 getParticleData(ib),
		 getParticleData(ic),
		 getParticleData(id)};
  // add to the Particle data pointer lists
  if(!pin[0]||!pin[1]||!pin[2]||!pin[3]) return;
  _particles.push_back(vector<PDPtr>(pin,pin+4));
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(!outgoing(_particles.back()[ix]->id())) {
      _outpart.push_back(_particles.back()[ix]);
      _ioutpart.push_back(_particles.back()[ix]->id());
    }
  }
  // add to the list of incoming particles  
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(_particles.back()[ix]->CC()) {
      if(!incoming(_particles.back()[ix]->CC()->id())) {
	_inpart.push_back(_particles.back()[ix]->CC());
	_iinpart.push_back(_particles.back()[ix]->CC()->id());
      }
    }
    else {
      if(!incoming(_particles.back()[ix]->id())) {
	_inpart.push_back(_particles.back()[ix]);
	_iinpart.push_back(_particles.back()[ix]->id());
      }
    }
  }
  // increment the size of the arrays
  ++_nsize;
}

// add particle to the list for a five point vertex
void VertexBase::add(int ia,int ib,int ic,int id, int ie) {
  if(_npoint!=5)
    throw HelicityConsistencyError() << "This is a " << _npoint 
				     << " vertex cannot add five particles" 
				     << Exception::abortnow;
  // add to the PDG code lists
  _iparticlea.push_back(ia);
  _iparticleb.push_back(ib);
  _iparticlec.push_back(ic);
  _iparticled.push_back(id);
  _iparticlee.push_back(ie);
  tPDPtr pin[5]={getParticleData(ia),
		 getParticleData(ib),
		 getParticleData(ic),
		 getParticleData(id),
		 getParticleData(ie)};
  // add to the Particle data pointer lists
  if(!pin[0]||!pin[1]||!pin[2]||!pin[3]||!pin[4]) return;
  _particles.push_back(vector<PDPtr>(pin,pin+5));
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(!outgoing(_particles.back()[ix]->id())) {
      _outpart.push_back(_particles.back()[ix]);
      _ioutpart.push_back(_particles.back()[ix]->id());
    }
  }
  // add to the list of incoming particles  
  for(unsigned int ix=0;ix<_particles.back().size();++ix) {
    if(_particles.back()[ix]->CC()) {
      if(!incoming(_particles.back()[ix]->CC()->id())) {
	_inpart.push_back(_particles.back()[ix]->CC());
	_iinpart.push_back(_particles.back()[ix]->CC()->id());
      }
    }
    else {
      if(!incoming(_particles.back()[ix]->id())) {
	_inpart.push_back(_particles.back()[ix]);
	_iinpart.push_back(_particles.back()[ix]->id());
      }
    }
  }
  // increment the size of the arrays
  ++_nsize;
}


// set the list of outgoing particles
void VertexBase::setOutgoing() {
  if(_outpart.size()!=0)
    throw HelicityConsistencyError() << "VertexBase::setOutgoing " 
				     << "Outgoing particles already set" 
				     << Exception::abortnow;
  for(unsigned int ix=0;ix<_particles.size();++ix) {
    for(unsigned int iy=0;iy<_particles[ix].size();++iy) {
      if(!outgoing(_particles[ix][iy]->id())) {
	_outpart.push_back(_particles[ix][iy]);
	_ioutpart.push_back(_particles[ix][iy]->id());
      }
    }
  } 
}

// set the list of incoming particles
void VertexBase::setIncoming() {
  if(_inpart.size()!=0) 
    throw HelicityConsistencyError() << "VertexBase::setIncoming " 
				     << "Outgoing particles already set" 
				     << Exception::abortnow;
  PDPtr temp;
  for(unsigned int ix=0;ix<_particles.size();++ix) {
    for(unsigned int iy=0;iy<_particles[ix].size();++iy) {
      temp =_particles[ix][iy]->CC();
      if(temp) {
	if(!incoming(temp->id())) {
	  _inpart.push_back(temp);
	  _iinpart.push_back(temp->id());
	}
      }
      else {
	if(!incoming(_particles[ix][iy]->id())) {
	  _inpart.push_back(_particles[ix][iy]);
	  _iinpart.push_back(_particles[ix][iy]->id());
	}
      }
    }
  }
}
 
// setup the lists of particles

// for a three point vertex
void VertexBase::setList(vector<int> ida ,vector<int> idb, vector<int> idc) {
  // check not already set
  if(_iparticlea.size()!=0)
    throw HelicityConsistencyError() << "VertexBase::setList the Particles have already"
				     << " been set" << Exception::abortnow;
  // check this is really a three point vertex
  else if(_npoint!=3)
    throw HelicityConsistencyError() << "VertexBase::setList not a three point vertex" 
				     << Exception::abortnow;
  // check all vectors have the same size
  else if(ida.size()!=ida.size()||ida.size()!=idc.size())
    throw HelicityConsistencyError() << "VertexBase::setList Particle lists have " 
				     << "different sizes" << Exception::abortnow;
  // set the arrays
  else {
    _iparticlea=ida; 
    _iparticleb=idb; 
    _iparticlec=idc;
  }
  _nsize = 0;
}

// for a four point vertex
void VertexBase::setList(vector<int> ida,vector<int> idb,
				vector<int> idc,vector<int> idd) {
  // check not already set
  if(_iparticlea.size()!=0)
    throw HelicityConsistencyError() << "VertexBase::setList the Particles have already"
				     << " been set" << Exception::abortnow;
  // check this is really a four point vertex
  else if(_npoint!=4)
    throw HelicityConsistencyError() << "VertexBase::setList not a four point vertex" 
				     << Exception::abortnow;
  // check all vectors have the same size
  else if(ida.size()!=ida.size()||ida.size()!=idc.size()||ida.size()!=idd.size())
    throw HelicityConsistencyError() << "VertexBase::setList Particle lists have " 
				     << "different sizes" << Exception::abortnow;
  // set the arrays
  else {
    _iparticlea=ida;
    _iparticleb=idb;
    _iparticlec=idc;
    _iparticled=idd;
  }
  _nsize = 0;
}

void VertexBase::setList(vector<int> ida,vector<int> idb, vector<int> idc,
			 vector<int> idd,vector<int> ide) {
  // check not already set
  if(_iparticlea.size()!=0)
    throw HelicityConsistencyError() << "VertexBase::setList the Particles have already"
				     << " been set" << Exception::abortnow;
  // check this is really a five point vertex
  else if(_npoint!=5)
    throw HelicityConsistencyError() << "VertexBase::setList not a five point vertex" 
				     << Exception::abortnow;
  // check all vectors have the same size
  else if(ida.size()!=ida.size()||ida.size()!=idc.size()||ida.size()!=idd.size()||
	  ida.size()!=ide.size())
    throw HelicityConsistencyError() << "VertexBase::setList Particle lists have " 
				     << "different sizes" << Exception::abortnow;
  // set the arrays
  else {
    _iparticlea=ida;
    _iparticleb=idb;
    _iparticlec=idc;
    _iparticled=idd;
    _iparticlee=ide;
  }
  _nsize = 0;
}
