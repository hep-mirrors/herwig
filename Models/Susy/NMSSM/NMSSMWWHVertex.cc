// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWWHVertex class.
//

#include "NMSSMWWHVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace Herwig::Helicity;

NMSSMWWHVertex::NMSSMWWHVertex() 
  : _couplast(0.), _q2last(), _mw(), _zfact(0.), _sw(0.), _sinb(0.),_cosb(0.) {
  int id[3]={25,35,45};
  vector<int> first,second,third;
  for(unsigned int ix=0;ix<3;++ix) {
    first .push_back(  24  );
    second.push_back( -24  );
    third .push_back(id[ix]);
    first .push_back(  23  );
    second.push_back(  23  );
    third .push_back(id[ix]);
  }
  setList(first,second,third);
}

void NMSSMWWHVertex::doinit() throw(InitException) {
  // SM parameters
  _theSM = generator()->standardModel();
  _mw=getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _sw = _theSM->sin2ThetaW();
  _zfact = 1./(1.-_sw);
  _sw=sqrt(_sw);
  // NMSSM parameters
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMWWHVertex::doinit()"
			  << Exception::runerror;
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWWHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb=sin(beta);
  _cosb=cos(beta);
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  // base class
  VVSVertex::doinit();
}

void NMSSMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _zfact << _sw << _sinb << _cosb << _mixS;
}

void NMSSMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _zfact >> _sw >> _sinb >> _cosb >> _mixS;
}

ClassDescription<NMSSMWWHVertex> NMSSMWWHVertex::initNMSSMWWHVertex;
// Definition of the static class description member.

void NMSSMWWHVertex::Init() {

  static ClassDocumentation<NMSSMWWHVertex> documentation
    ("The NMSSMWWHVertex class implements the coupling of two electroweak gauge"
     " bosons with the Higgs bosons of the NMSSM");

}

void NMSSMWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr c) {
  // gauge bosons
  int ibos=abs(a->id());
  // higgs
  int ihiggs = (c->id()-25)/10;
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = sqrt(4.0*Constants::pi*alpha)*_mw/_sw*UnitRemoval::InvE;
    _q2last=q2;
  }
  // higgs mixing factor
  Complex hmix = _cosb*(*_mixS)(ihiggs,0)+_sinb*(*_mixS)(ihiggs,1);
  // coupling
  if(ibos==24)      setNorm(_couplast*hmix);
  else if(ibos==23) setNorm(_couplast*hmix*_zfact);
  else {
    throw HelicityConsistencyError() << "SMWWHVertex::setCoupling "
				     << "Invalid particles in WWH Vertex" 
				     << Exception::warning;
    setNorm(0.);
  }
}
