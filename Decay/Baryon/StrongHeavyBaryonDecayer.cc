// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StrongHeavyBaryonDecayer class.
//

#include "StrongHeavyBaryonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "StrongHeavyBaryonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

StrongHeavyBaryonDecayer::~StrongHeavyBaryonDecayer() {}

void StrongHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _gsigma_clambda_cpi << _gxistar_cxi_cpi << _flambda_c1sigma_cpi 
     << _flambda_c1starsigma_cpi  << _gsigma_blambda_bpi << _gxistar_bxi_bpi 
     << _flambda_b1sigma_bpi << _flambda_b1starsigma_bpi <<  _incoming << _outgoingB 
     << _outgoingM << _maxweight << _A1 << _A2 << _A3 << _B1 << _B2 << _B3;
}

void StrongHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _gsigma_clambda_cpi >> _gxistar_cxi_cpi >> _flambda_c1sigma_cpi 
     >> _flambda_c1starsigma_cpi  >> _gsigma_blambda_bpi >> _gxistar_bxi_bpi 
     >> _flambda_b1sigma_bpi >> _flambda_b1starsigma_bpi >>  _incoming >> _outgoingB 
     >> _outgoingM >> _maxweight >> _A1 >> _A2 >> _A3 >> _B1 >> _B2 >> _B3;
}

bool StrongHeavyBaryonDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed=false;
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;
  do
    {
      if(id0==_incoming[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incoming[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){allowed=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&!allowed);
  return allowed;
}

ParticleVector StrongHeavyBaryonDecayer::decay(const DecayMode & dm,
						 const Particle & parent) const {
  int imode=-1;
  int id=parent.id();
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  unsigned int ix=0;bool cc;
  do 
    {
      if(id==_incoming[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incoming[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;cc=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_incoming.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}

ClassDescription<StrongHeavyBaryonDecayer> StrongHeavyBaryonDecayer::initStrongHeavyBaryonDecayer;
// Definition of the static class description member.

void StrongHeavyBaryonDecayer::Init() {

  static ClassDocumentation<StrongHeavyBaryonDecayer> documentation
    ("There is no documentation for the \\classname{StrongHeavyBaryonDecayer} class");

}

// couplings for spin-1/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
halfHalfScalarCoupling(int imode,Complex & A,Complex & B) const
{
  A =_A1[imode];B=_B1[imode];  
}

// couplings for spin-1/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
halfThreeHalfScalarCoupling(int imode,Complex& A,Complex& B) const
{
  A =_A1[imode];B=_B1[imode];  
}


}
