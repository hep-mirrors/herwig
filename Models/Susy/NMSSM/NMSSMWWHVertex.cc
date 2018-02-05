// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWWHVertex class.
//

#include "NMSSMWWHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMWWHVertex::NMSSMWWHVertex() 
  : _couplast(0.), _q2last(), _mw(), _zfact(0.), _sinb(0.),_cosb(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMWWHVertex::doinit() {
  int id[3]={25,35,45};
  // PDG codes for the particles in the vertex
  for(unsigned int ix=0;ix<3;++ix) {
    // Higgs WW
    addToList( 24, -24, id[ix] );
    //Higgs ZZ
    addToList( 23, 23, id[ix] );
  }
  // SM parameters
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  double sw2 = sin2ThetaW();
  _zfact = 1./(1.-sw2);
  // NMSSM parameters
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
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
  // base class
  VVSVertex::doinit();
}

void NMSSMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw,GeV) << _zfact << _sinb << _cosb << _mixS;
}

void NMSSMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw,GeV) >> _zfact >> _sinb >> _cosb >> _mixS;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMWWHVertex,VVSVertex>
describeHerwigNMSSMWWHVertex("Herwig::NMSSMWWHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMWWHVertex::Init() {

  static ClassDocumentation<NMSSMWWHVertex> documentation
    ("The NMSSMWWHVertex class implements the coupling of two electroweak gauge"
     " bosons with the Higgs bosons of the NMSSM");

}
//calulate couplings
void NMSSMWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr c) {
  // ID of gauge bosons
  int ibos=abs(a->id());
  // ID of Higgs
  int ihiggs = (c->id()-25)/10;
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2)*_mw*UnitRemoval::InvE;
    _q2last=q2;
  }
  // higgs mixing factor
  Complex hmix = _cosb*(*_mixS)(ihiggs,0)+_sinb*(*_mixS)(ihiggs,1);
  // couplings
  if(ibos==24)      norm(_couplast*hmix);
  else if(ibos==23) norm(_couplast*hmix*_zfact);
  else assert(false);
}
