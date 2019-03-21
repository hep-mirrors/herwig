// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGGHVertex class.
//

#include "NMSSMGGHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/Susy/NMSSM/NMSSM.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

NMSSMGGHVertex::NMSSMGGHVertex() : _sw(0.), _cw(0.), _mw(0.*MeV),
	_mz(0.*MeV),_lambdaVEV(0.*MeV), _lambda(0.), _v1(0.*MeV),
	_v2(0.*MeV), _triTp(0.*MeV), _triBt(0.*MeV),
	_sb(0.), _cb(0.), _masslast(make_pair(0.*MeV,0.*MeV)),
	_q2last(0.*MeV2), _couplast(0.), _coup(0.),
    _hlast(0), _recalc(true) {
  orderInGem(1);
  orderInGs(2);
}

void NMSSMGGHVertex::doinit()  {
  addToList(21,21,25);
  addToList(21,21,35);
  addToList(21,21,36);
  addToList(21,21,45);
  addToList(21,21,46);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMGGHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1. - sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _top = getParticleData(6);
  _bt = getParticleData(5);
  
  //NMSSM parameters
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixQt = nmssm->stopMix();
  _mixQb = nmssm->sbottomMix();
  
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  
  _v1 = sqrt(2.)*_mw*_cb;
  _v2 = sqrt(2.)*_mw*_sb;
 
  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();

  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();

  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(6, Energy());
  masses[0] = getParticleData(6)->mass();
  masses[1] = getParticleData(5)->mass();

  masses[2] = getParticleData(1000005)->mass();
  masses[3] = getParticleData(2000005)->mass();

  masses[4] = getParticleData(1000006)->mass();
  masses[5] = getParticleData(2000006)->mass();

  type.resize(6, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1Half;
  couplings.resize(6);

  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void NMSSMGGHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
     << ounit(_lambdaVEV,GeV) << _lambda << ounit(_v1,GeV) << ounit(_v2,GeV)
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) 
     << _top << _bt << _mixS << _mixP << _mixQt << _mixQb << _sb << _cb; 
}


void NMSSMGGHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> iunit(_lambdaVEV,GeV) >> _lambda >> iunit(_v1,GeV) >> iunit(_v2,GeV)
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) 
     >> _top >> _bt >> _mixS >> _mixP >> _mixQt >> _mixQb >> _sb >> _cb; 
}

ClassDescription<NMSSMGGHVertex> NMSSMGGHVertex::initNMSSMGGHVertex;
// Definition of the static class description member.

void NMSSMGGHVertex::Init() {

  static ClassDocumentation<NMSSMGGHVertex> documentation
    ("The effective coupling of a higgs to a pair of gluons in the "
     "NMSSM.");

}

void NMSSMGGHVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
				 tcPDPtr p3) {			 
  long hid(p3->id());
  if( q2 != _q2last ) {
    Looptools::clearcache();
    _couplast = sqr(strongCoupling(q2));
    _coup = weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  norm(_couplast*_coup);
  // scalar higgs bosons  
  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    if( hid % 5 == 0 ) {
      // location of the higgs
      int iloc = (hid - 25)/10;
      // 6 particles in the loop
      setNParticles(6);
      // top and bottom quark masses
      Energy mt = _theSM->mass(q2, _top);
      Energy mb = _theSM->mass(q2,  _bt);
      Complex c(0.);
      // couplings for the top quark loop
      c = -0.25*mt*(*_mixS)(iloc, 1)/_sb/_mw;
	
      couplings[0] = make_pair(c,c);
      masses[0] = mt;
      // couplings for the bottom quark loop
      c = -0.25*mb*(*_mixS)(iloc, 0)/_cb/_mw;
	
      couplings[1] = make_pair(c,c);	
      masses[1] = mb;
      // sbottoms
      double f1 = mb/_mw/_cb;
      complex<Energy>  f2 = 0.5*_mz/_cw*
	( - _cb*(*_mixS)(iloc,0) + _sb*(*_mixS)(iloc,1));
      complex<Energy> cpl;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw)/3.)*(*_mixQb)(ix, 0)*(*_mixQb)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)/3.)
	  - f1*mb*(*_mixS)(iloc,0)
	  *((*_mixQb)(ix, 0)*(*_mixQb)(ix, 0) + (*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_v2*(*_mixS)(iloc,2)/_coup 
		    +  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(ix, 1)*(*_mixQb)(ix, 0)
						 + (*_mixQb)(ix, 0)*(*_mixQb)(ix, 1));

	couplings[2+ix] = make_pair(Complex(0.5*cpl*UnitRemoval::InvE),
				    Complex(0.5*cpl*UnitRemoval::InvE)); 
      }
      // stop
      f1 = mt/_mw/_sb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl   =+f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
		     + 4.*sqr(_sw)*(*_mixQt)(ix, 1)*(*_mixQt)(ix, 1)/3.)
	  - f1*mt*(*_mixS)(iloc,1)
	  *((*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
	    + (*_mixQt)(ix, 1)*(*_mixQt)(ix, 1))  
	  -  0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,0) - _lambda*_v1*(*_mixS)(iloc,2)/_coup
		     + _triTp*(*_mixS)(iloc,1))*((*_mixQt)(ix, 1)*(*_mixQt)(ix, 0)
						 + (*_mixQt)(ix, 0)*(*_mixQt)(ix, 1));

	couplings[4+ix] = make_pair(Complex(0.5*cpl*UnitRemoval::InvE),
				    Complex(0.5*cpl*UnitRemoval::InvE));
      }
    }
    // pseudoscalar higgs bosons	
    else {
      // location of the higgs
      int iloc = (hid - 36)/10;
      // 2 particles in the loop
      setNParticles(2);
      // top and bottom quark masses
      Energy mt = _theSM->mass(q2, _top);
      Energy mb = _theSM->mass(q2,  _bt);
      Complex c(0.);
      // top quark couplings
      c = Complex(0.,-1.)*0.25*mt*(*_mixP)(iloc, 1)/_sb/_mw;
	  	  
      couplings[0] = make_pair(-c,c);
      masses[0] = mt;
      // bottom quark couplings
      c = Complex(0., -1.)*0.25*mb*(*_mixP)(iloc, 0)/_cb/_mw;
	  
      couplings[1] = make_pair(-c,c);	
      masses[1] = mb;
    }
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
