// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VertexBase class.
//

#include "VertexBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

namespace Herwig{
namespace Helicity{
using namespace ThePEG;
using namespace ThePEG::Helicity;
    
VertexBase::~VertexBase() {}
    
void VertexBase::persistentOutput(PersistentOStream & os) const {
  os << _npoint << _nsize << _ispin << _inpart << _iinpart << _outpart << _ioutpart 
     << _iparticlea << _iparticleb << _iparticlec << _iparticled << _iparticlee
     << _particlea << _particleb << _particlec << _particled << _particlee;
}

void VertexBase::persistentInput(PersistentIStream & is, int) {
  is >> _npoint >> _nsize >> _ispin >> _inpart >> _iinpart >> _outpart >> _ioutpart 
     >> _iparticlea >> _iparticleb >> _iparticlec >> _iparticled >> _iparticlee
     >> _particlea >> _particleb >> _particlec >> _particled >> _particlee;
}
    
ClassDescription<VertexBase> VertexBase::initVertexBase;
// Definition of the static class description member.
  
void VertexBase::Init() {
  
  static Parameter<VertexBase,int> interfacenpoint
    ("NPoint",
     "The number of extermal particles interacting at the Vertex.",
     &VertexBase::_npoint, 3, 3, 5, false, false, true);
  
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
    ("The \\classname{VertexBase} class is designed to be the base class"
     "of all vertices in Herwig++");
  
}

// find particles with a given id    
vector<PDPtr> VertexBase::search(int iloc,int idd)
{
  vector<PDPtr> out;
  if(iloc<0||iloc>=_npoint)
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
	  if(_npoint<=3)
	    {
	      out.push_back(_particlea[ix]);
	      out.push_back(_particleb[ix]);
	      out.push_back(_particlec[ix]);
	    }
	  if(_npoint<=4){out.push_back(_particled[ix]);}
	  if(_npoint==5){out.push_back(_particlee[ix]);}
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
   
}
}

