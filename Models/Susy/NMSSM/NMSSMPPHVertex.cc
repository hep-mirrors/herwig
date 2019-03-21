// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMPPHVertex class.
//

#include "NMSSMPPHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Models/Susy/NMSSM/NMSSM.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

NMSSMPPHVertex::NMSSMPPHVertex() 
  : _sw(0.), _cw(0.), _mw(0.*MeV),
    _mz(0.*MeV),_lambdaVEV(0.*MeV), _lambda(0.),
    _triTp(0.*MeV), _triBt(0.*MeV),
    _sb(0.), _cb(0.), 
    _kappa(0.),_vu(ZERO),_vd(ZERO),_s(ZERO),_theAl(ZERO),
    _masslast(make_pair(0.*MeV,0.*MeV)),_q2last(0.*MeV2),
    _couplast(0.), _coup(0.), _hlast(0), _recalc(true) {
  orderInGem(3);
  orderInGs(0);
}

void NMSSMPPHVertex::doinit()  {
  addToList(22,22,25);
  addToList(22,22,35);
  addToList(22,22,36);
  addToList(22,22,45);
  addToList(22,22,46);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMPPHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1. - sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _top = getParticleData(6);
  _bt  = getParticleData(5);
  _tau = getParticleData(15);
  
  //NMSSM parameters
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixQt = nmssm->stopMix();
  _mixQb = nmssm->sbottomMix();
  _mixLt = nmssm->stauMix();
  
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);

  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();
  
  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();
  _triTa = nmssm->tauTrilinear();

  _vd = sqrt(2)*_mw*_cb;
  _vu = sqrt(2)*_mw*_sb;
  _s  = _lambdaVEV/_lambda;
  _theAl = nmssm->trilinearLambda();
  _kappa = nmssm->kappa();

  _mixU = nmssm->charginoUMix();
  _mixV = nmssm->charginoVMix();

  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(13, Energy());
  masses[ 0] = getParticleData( 6)->mass();
  masses[ 1] = getParticleData( 5)->mass();
  masses[ 2] = getParticleData(15)->mass();
  masses[ 3] = getParticleData(ParticleID::SUSY_chi_1plus)->mass();
  masses[ 4] = getParticleData(ParticleID::SUSY_chi_2plus)->mass();
  masses[ 5] = _mw;
  masses[ 6] = getParticleData(ParticleID::Hplus)->mass();
  masses[ 7] = getParticleData(1000005)->mass();
  masses[ 8] = getParticleData(2000005)->mass();
  masses[ 9] = getParticleData(1000006)->mass();
  masses[10] = getParticleData(2000006)->mass();
  masses[11] = getParticleData(1000015)->mass();
  masses[12] = getParticleData(2000015)->mass();
  type.resize(13, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1Half;
  type[2] = PDT::Spin1Half;
  type[3] = PDT::Spin1Half;
  type[4] = PDT::Spin1Half;
  type[5] = PDT::Spin1;
  couplings.resize(13);
  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void NMSSMPPHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
     << ounit(_lambdaVEV,GeV) << _lambda
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) << ounit(_triTa,GeV)
     << _top << _bt << _tau << _mixS << _mixP << _mixU << _mixV
     << _mixQt << _mixQb << _mixLt << _sb << _cb <<  _kappa
     << ounit(_vu,GeV) << ounit(_vd,GeV) << ounit(_s,GeV) << ounit(_theAl,GeV);
}

void NMSSMPPHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> iunit(_lambdaVEV,GeV) >> _lambda
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) >> iunit(_triTa,GeV)
     >> _top >> _bt >> _tau >> _mixS >> _mixP >> _mixU >> _mixV
     >> _mixQt >> _mixQb >> _mixLt >> _sb >> _cb >> _kappa 
     >> iunit(_vu,GeV) >> iunit(_vd,GeV) >> iunit(_s,GeV) >> iunit(_theAl,GeV);
}

ClassDescription<NMSSMPPHVertex> NMSSMPPHVertex::initNMSSMPPHVertex;
// Definition of the static class description member.

void NMSSMPPHVertex::Init() {

  static ClassDocumentation<NMSSMPPHVertex> documentation
    ("The effective coupling of a higgs to a pair of gluons in the "
     "NMSSM.");

}

void NMSSMPPHVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
				 tcPDPtr p3) {			 
  long hid(p3->id());
  double rt = sqrt(0.5);
  if( q2 != _q2last ) {
    Looptools::clearcache();
    _couplast = sqr(electroMagneticCoupling(q2));
    _coup = weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  norm(_couplast*_coup);
  // scalar higgs bosons  
  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    // top and bottom quark masses
    Energy mt   = _theSM->mass(q2, _top);
    Energy mb   = _theSM->mass(q2,  _bt);
    Energy mtau = _theSM->mass(q2, _tau);
    // scalar
    if( hid % 5 == 0 ) {
      // location of the higgs
      int iloc = (hid - 25)/10;
      // 6 particles in the loop
      setNParticles(13);
      Complex c(0.);
      // couplings for the top quark loop
      c = -1.5*sqr(_theSM->eu())*  mt*(*_mixS)(iloc, 1)/_sb/_mw;
      couplings[0] = make_pair(c,c);
      masses[0] = mt;
      // couplings for the bottom quark loop
      c = -1.5*sqr(_theSM->ed())*  mb*(*_mixS)(iloc, 0)/_cb/_mw;
      couplings[1] = make_pair(c,c);	
      masses[1] = mb;
      // couplings for the tau lepton loop
      c = -0.5*sqr(_theSM->ee())*mtau*(*_mixS)(iloc, 0)/_cb/_mw;
      couplings[2] = make_pair(c,c);	
      masses[2] = mtau;
      // charginos
      for(unsigned int ic=0;ic<2;++ic) {
	c = -_lambda/_coup*rt*(*_mixS)(iloc,2)*(*_mixU)(ic,1)*(*_mixV)(ic,1)
	  -rt*((*_mixS)(iloc,0)*(*_mixU)(ic,1)*(*_mixV)(ic,0) +
	       (*_mixS)(iloc,1)*(*_mixU)(ic,0)*(*_mixV)(ic,1));
	couplings[3+ic] = make_pair(c,c);
      }
      // W boson
      c = Complex(UnitRemoval::InvE*_mw*
		  (_cb*(*_mixS)(iloc,0)+_sb*(*_mixS)(iloc,1)));
      couplings[5] = make_pair(c,c);
      // charged Higgs
      complex<Energy> cpl;
      cpl = sqr(_lambda)*rt*2.*(_s*((*_mixS)(iloc,2)*sqr(_cb) + (*_mixS)(iloc,2)*sqr(_sb))
				- (_vu*(*_mixS)(iloc,0)/_coup + 
				   _vd*(*_mixS)(iloc,1)/_coup)*_sb*_cb)
	+_lambda*_sb*_cb*2*(*_mixS)(iloc,2)*(_kappa*_s/rt + rt*_theAl)
	+ sqr(_coup)*0.5*rt*sqr(_sw)/sqr(_cw)*((_vu*(*_mixS)(iloc,1)/_coup - 
						_vd*(*_mixS)(iloc,0)/_coup)*sqr(_cb) + 
					       (_vd*(*_mixS)(iloc,0)/_coup -
						_vu*(*_mixS)(iloc,1)/_coup)*sqr(_sb))
	+ sqr(_coup)*0.5*rt*(_vu*((*_mixS)(iloc,1)*sqr(_cb) + 
				  (*_mixS)(iloc,1)*sqr(_sb) + 
				  2.*(*_mixS)(iloc,0)*_cb*_sb)/_coup 
			     + _vd*((*_mixS)(iloc,0)*sqr(_cb) + 
				    (*_mixS)(iloc,0)*sqr(_sb) 
				    + 2.*(*_mixS)(iloc,1)*_sb*_cb)/_coup);
      cpl /= -_coup;
      couplings[6] = make_pair(Complex(cpl*UnitRemoval::InvE),
			       Complex(cpl*UnitRemoval::InvE));
      // sbottoms
      double f1 = mb/_mw/_cb;
      complex<Energy>  f2 = 0.5*_mz/_cw*
	( - _cb*(*_mixS)(iloc,0) + _sb*(*_mixS)(iloc,1));
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw)/3.)*(*_mixQb)(ix, 0)*(*_mixQb)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)/3.)
	  - f1*mb*(*_mixS)(iloc,0)
	  *((*_mixQb)(ix, 0)*(*_mixQb)(ix, 0) + (*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_vu*(*_mixS)(iloc,2)/_coup 
		    +  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(ix, 1)*(*_mixQb)(ix, 0)
						 + (*_mixQb)(ix, 0)*(*_mixQb)(ix, 1));
	cpl *= 3.*sqr(_theSM->ed());
	couplings[7+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),Complex(cpl*UnitRemoval::InvE)); 
      }
      // stop
      f1 = mt/_mw/_sb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl   =+f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
		     + 4.*sqr(_sw)*(*_mixQt)(ix, 1)*(*_mixQt)(ix, 1)/3.)
	  - f1*mt*(*_mixS)(iloc,1)
	  *((*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
	    + (*_mixQt)(ix, 1)*(*_mixQt)(ix, 1))  
	  -  0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,0) - _lambda*_vd*(*_mixS)(iloc,2)/_coup
		     + _triTp*(*_mixS)(iloc,1))*((*_mixQt)(ix, 1)*(*_mixQt)(ix, 0)
						 + (*_mixQt)(ix, 0)*(*_mixQt)(ix, 1));
	cpl *= 3.*sqr(_theSM->eu());
	couplings[9+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),
				    Complex(cpl*UnitRemoval::InvE));
      } // sbottoms
      f1 = mtau/_mw/_cb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw))*(*_mixLt)(ix, 0)*(*_mixLt)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixLt)(ix, 1)*(*_mixLt)(ix, 1))
	  - f1*mtau*(*_mixS)(iloc,0)
	  *((*_mixLt)(ix, 0)*(*_mixLt)(ix, 0) + (*_mixLt)(ix, 1)*(*_mixLt)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_vu*(*_mixS)(iloc,2)/_coup 
		    +  _triTa*(*_mixS)(iloc,0))*((*_mixLt)(ix, 1)*(*_mixLt)(ix, 0)
						 + (*_mixLt)(ix, 0)*(*_mixLt)(ix, 1));
	cpl *= sqr(_theSM->ee());
	couplings[11+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),
				     Complex(cpl*UnitRemoval::InvE)); 
      }
    }
    // pseudoscalar higgs bosons	
    else {
      // location of the higgs
      int iloc = (hid - 36)/10;
      // 2 particles in the loop
      setNParticles(5);
      Complex c(0.);
      // top quark couplings
      c = Complex(0., 1.)*1.5*sqr(_theSM->eu())*  mt*(*_mixP)(iloc, 1)/_sb/_mw;
      couplings[0] = make_pair(c,-c);
      masses[0] = mt;
      // bottom quark couplings
      c = Complex(0., 1.)*1.5*sqr(_theSM->ed())*  mb*(*_mixP)(iloc, 0)/_cb/_mw;
      couplings[1] = make_pair(c,-c);	
      masses[1] = mb;
      // tau lepton couplings
      c = Complex(0., 1.)*0.5*sqr(_theSM->ee())*mtau*(*_mixP)(iloc, 0)/_cb/_mw;
      couplings[2] = make_pair(c,-c);	
      masses[2] = mtau;
      // charginos
      for(unsigned int ic=0;ic<2;++ic) {
	c = Complex(0,-1.0)*
	  (_lambda/_coup*rt*(*_mixP)(iloc,2)*(*_mixU)(ic,1)*(*_mixV)(ic,1)
	   -rt*((*_mixP)(iloc,0)*(*_mixU)(ic,1)*(*_mixV)(ic,0)
		     + (*_mixP)(iloc,1)*(*_mixU)(ic,0)*(*_mixV)(ic,1)));
			  couplings[3+ic] = make_pair(-c,c); 

      }
    }
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
