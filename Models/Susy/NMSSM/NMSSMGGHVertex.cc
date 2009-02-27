// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGGHVertex class.
//

#include "NMSSMGGHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/Susy/NMSSM/NMSSM.h"

using namespace Herwig;

NMSSMGGHVertex::NMSSMGGHVertex() : 
  _sw(0.), _cw(0.), _mw(), _mz(), _sb(0.), _cb(0.), 
  _masslast(make_pair(0.*MeV,0.*MeV)), _q2last(0.*MeV2), _couplast(0.), 
  _hlast(0), _recalc(true) {
  //PDG codes for particles at vertices
  vector<long> gluon(5,21), third(5);
  third[0] = 25;
  third[1] = 35;
  third[2] = 36;
  third[3] = 45;
  third[4] = 46;
  setList(gluon, gluon, third);
}

void NMSSMGGHVertex::doinit() {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMGGHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(_theSM->sin2ThetaW());
  _cw = sqrt(1. - _theSM->sin2ThetaW());
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

  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(5, Energy());
  masses[1] = getParticleData(1000005)->mass();
  masses[2] = getParticleData(2000005)->mass();
  masses[3] = getParticleData(1000006)->mass();
  masses[4] = getParticleData(2000006)->mass();
  type.resize(5, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  couplings.resize(5);

  orderInGem(1);
  orderInGs(2);

  VVSLoopVertex::doinit();
}

void NMSSMGGHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
     << _top << _bt << _mixS << _mixP << _mixQt << _mixQb << _sb 
     << _cb; 
}

void NMSSMGGHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> _top >> _bt >> _mixS >> _mixP >> _mixQt >> _mixQb >> _sb
     >> _cb; 
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
    _masslast.first = _theSM->mass(q2, _bt);
    _masslast.second = _theSM->mass(q2, _top);
    //factor of 0.5 factored from each EW vertex
    _couplast = 0.5*sqr(strongCoupling(q2))*weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  setNorm(_couplast);

  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    if( hid % 5 == 0 ) {
      setNParticles(5);
      int iloc = (hid - 25)/10;
      //top quark
      masses[0] = _masslast.second;
      Complex cpl = -_masslast.second*(*_mixS)(iloc, 1)/_sb/_mw;
      couplings[0].first = cpl;
      couplings[0].second = cpl;
      //~b_1
      double f1 = _masslast.first/_mw/_cb;
      complex<Energy> f2 = 0.5*_mz*(_cb*(*_mixS)(iloc,0) 
				    - _sb*(*_mixS)(iloc,1))/_cw;
      Complex mix = (*_mixQb)(0, 0);
      cpl = ( f2*(1. - 2.*sqr(_sw)/3.) 
	      - f1*_masslast.first*(*_mixS)(iloc,0) )*mix*mix 
	* UnitRemoval::InvE;
      couplings[1].first = cpl; couplings[1].second = cpl; 
      //~b_2
      mix = (*_mixQb)(1, 1);
      cpl = -( f2*2.*sqr(_sw)/3. + f1*_masslast.first*(*_mixS)(iloc,0) )*mix*mix
	* UnitRemoval::InvE;
      couplings[2].first = cpl; couplings[2].second = cpl; 
      //~t_1
      f1 = _masslast.second/_mw/_sb;
      mix = (*_mixQt)(0, 0);
      cpl = -( f2*(1. - 4.*sqr(_sw)/3.) + f1*_masslast.second*(*_mixS)(iloc,1) )
	*mix*mix*UnitRemoval::InvE;
      couplings[3].first = cpl; couplings[3].second = cpl; 
      //~t_2
      mix = (*_mixQt)(1, 1);
      cpl = -( f2*4.*sqr(_sw)/3. + f1*_masslast.second*(*_mixS)(iloc,1) )
	*mix*mix*UnitRemoval::InvE;
      couplings[4].first = cpl; couplings[4].second = cpl; 
    }
    else {
      setNParticles(1);
      int iloc = (hid - 36)/10;
      //top quark
      masses[0] = _masslast.second;
      Complex cpl = Complex(0., 1.)*_masslast.second*(*_mixP)(iloc, 1)/_sb/_mw;
      couplings[0].first = cpl;
      couplings[0].second = -cpl;
    }
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
