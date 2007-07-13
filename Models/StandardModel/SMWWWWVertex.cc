// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWWVertex class.
//

#include "SMWWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void SMWWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _gamma << _Z0 << _wplus << _wminus
     << _vfact  << _sw2 << _cw2;
}

void SMWWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _gamma >> _Z0 >> _wplus >> _wminus
     >> _vfact >> _sw2 >> _cw2;
  _couplast = 0.;
  _q2last = 0.*GeV2;
}

ClassDescription<SMWWWWVertex>SMWWWWVertex::initSMWWWWVertex;
// Definition of the static class description member.

void SMWWWWVertex::Init() {
  static ClassDocumentation<SMWWWWVertex> documentation
    ("The SMWWWWVertex class is the implementation of the"
     " Standard Model quartic electroweka gauge boson coupling.");
  
}


// couplings for the WWWW vertex
void SMWWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				      tcPDPtr c,tcPDPtr d)
{
  // id's of the particles
  int id[4]={a->id(),b->id(),c->id(),d->id()};
  // order the particles
  int ngamma(0),nz(0);
  int iorder[4];
  for(int ix=0;ix<4;++ix)
    {
      if (id[ix]==22)
	ngamma+=1;
      else if (id[ix]==23)
	nz+=1;
    }
  // if photons or Z's
  if(ngamma!=0 || nz!=0)
    {
      int iy=0;
      // put the photons first
      for(int ix=0;iy<ngamma&&ix<4;++ix) {
	if(id[ix]==22) {
	  iorder[iy]=ix;
	  ++iy;
	}
      }
      // then the Z bosons
      for(int ix=0;iy<ngamma+nz&&ix<4;++ix) {
	if(id[ix]==23) {
	  iorder[iy]=ix;
	  ++iy;
	}
      }
      // then the W+
      for(int ix=0;iy<3&&ix<4;++ix) {
	if(id[ix]==24) {
	  iorder[iy]=ix;
	  ++iy;
	}
      }
      if (iy!=3)
        {
	  throw HelicityConsistencyError() << "SMWWWWVertex::setCoupling"
					   << " Error setting order" 
					   << Exception::warning;
	  setNorm(0.);
	  return;
	}
      // finally the W-
      for(int ix=0;iy<4&&ix<4;++ix) {
	if(id[ix]==-24) {
	  iorder[iy]=ix;
	  ++iy;
	}
      }
      if(iy!=4)
        {
	  throw HelicityConsistencyError() << "SMWWWWVertex::setCoupling"
					   << " Error setting order" 
					   << Exception::warning;
	  setNorm(0.);return;
	}
    }
  else
    {
      int iy=0;
      // first the W+
      for(int ix=0;iy<3&&ix<4;++ix) {
	if(id[ix]==24) {
	  iorder[iy]=ix;++iy;
	}
      }
      if(iy!=2)
        {
	  throw HelicityConsistencyError() << "SMWWWWVertex::setCoupling"
					   << " Error setting order" 
					   << Exception::warning;
	  setNorm(0.);return;
	}
      // finally the W-
      for(int ix=0;iy<4&&ix<4;++ix) {
	if(id[ix]==-24) {
	  iorder[iy]=ix;
	  ++iy;
	}
      }
      if(iy!=4)
        {
	  throw HelicityConsistencyError() << "SMWWWWVertex::setCoupling"
					   << " Error setting order" 
					   << Exception::warning;
	  setNorm(0.);return;
	}
      setIntermediate(_gamma,_Z0,_sw2,_cw2);
    }
  setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);
  setType(2);
  // first the overall normalisation
  if(q2!=_q2last)
    {
      double alpha = _theSM->alphaEM(q2);
      _couplast = 4.0*Constants::pi*alpha;
      _q2last=q2;
    }
  // id's of the first two particles
  int ida(0),idb(0);
  if(iorder[0]==0){ida = abs(a->id());}
  else if(iorder[0]==1){ida = abs(b->id());}
  else if(iorder[0]==2){ida = abs(c->id());}
  else if(iorder[0]==3){ida = abs(d->id());}
  if(iorder[1]==0){idb = abs(a->id());}
  else if(iorder[1]==1){idb = abs(b->id());}
  else if(iorder[1]==2){idb = abs(c->id());}
  else if(iorder[1]==3){idb = abs(d->id());}
  // WWWW coupling
  if(ida==24){setNorm(_vfact[0]*_couplast);}
  // ZZWW coupling
  else if(ida==23&&idb==23){setNorm(_vfact[1]*_couplast);}
  // gamma gamma WW coupling
  else if(ida==22&&idb==22){setNorm(_couplast);}
  // gamma  Z WW coupling
  else if(ida==22&&idb==23){setNorm(_vfact[3]*_couplast);}
  else
    {
      throw HelicityConsistencyError() << "SMWWWWVertex::setCoupling"
				       << "unknown particles" 
				       << Exception::warning;
      setNorm(0.);return;
    }
}
 
}
