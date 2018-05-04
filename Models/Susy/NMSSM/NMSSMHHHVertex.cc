// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHHHVertex class.
//

#include "NMSSMHHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMHHHVertex::NMSSMHHHVertex() : _mw(0.*MeV), _mz(0.*MeV), _sw2(0.),
				   _cw(0.), _lambda(0.), _kappa(0.) ,
				   _lambdaVEV(0.*MeV), _theAl(0.*MeV),
				   _theAk(0.*MeV), _sb(0.), _cb(0.),
				   _s2b(0.), _c2b(0.), _vu(0.*MeV),
				   _vd(0.*MeV), _s(0.*MeV), _q2last(0.*MeV2),
				   _glast(0.), _MQ3(0.*MeV), _MU2(0.*MeV),
				   _includeRadiative(false) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMHHHVertex::doinit() {
  // PDG codes for the particles in vertex _vd
  //CP-even Higgs
  addToList(25, 35, 45);
  for( unsigned int i = 25; i <= 45; i += 10 ) {
    addToList(i, i, 25);
    addToList(i, i, 35);
    addToList(i, i, 45);
    //Charged Higgs
    addToList(i, 37, -37);
    //CP-odd Higgs
    addToList(i, 36, 36);
    addToList(i, 36, 46);
    addToList(i, 46, 36);
    addToList(i, 46, 46);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if( !nmssm ) 
    throw InitException() << "NMSSMHHHVertex::doinit - The model object is"
			  << "not the NMSSM object."
			  << Exception::runerror;  
  //SM parameters
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _sw2 = sin2ThetaW(); 
  _cw = sqrt(1. - _sw2);
  //NMSSM parameters
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  if( !_mixS || !_mixP ) 
    throw InitException() << "NMSSMHHHVertex::doinit - One of the mixing matrix "
			  << "pointers is null, cannot continue. S: "
			  << _mixS << "  P: " << _mixP << Exception::runerror;
  _lambda = nmssm->lambda();
  _kappa = nmssm->kappa();
  _lambdaVEV = nmssm->lambdaVEV();
  _theAl = nmssm->trilinearLambda();
  _theAk = nmssm->trilinearKappa();
  _MQ3 = nmssm->MQ3();  
  _MU2 = nmssm->MU2(); 
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  _vd = sqrt(2)*_mw*_cb;
  _vu = sqrt(2)*_mw*_sb;
  _s  = _lambdaVEV/_lambda;
  SSSVertex::doinit();
}

void NMSSMHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw, GeV) << ounit(_mz,GeV) 
     << _sw2 << _cw <<  _lambda << _includeRadiative
     << _kappa <<  ounit(_lambdaVEV,GeV) <<  ounit(_theAl, GeV) 
     << ounit(_theAk,GeV) <<  _sb <<  _cb << _s2b <<  _c2b
     << ounit(_vu,GeV) << ounit(_vd,GeV) << ounit(_s,GeV) << _mixS << _mixP
     << ounit(_MQ3,GeV) << ounit(_MU2,GeV) << _theSM;
}

void NMSSMHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw, GeV) >> iunit(_mz,GeV) 
     >> _sw2 >> _cw >>  _lambda >> _includeRadiative
     >> _kappa >>  iunit(_lambdaVEV,GeV) >>  iunit(_theAl, GeV) 
     >> iunit(_theAk,GeV) >>  _sb >>  _cb >> _s2b >>  _c2b
     >> iunit(_vu,GeV) >> iunit(_vd,GeV) >> iunit(_s,GeV)>> _mixS >> _mixP
     >> iunit(_MQ3,GeV) >> iunit(_MU2,GeV)  >> _theSM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMHHHVertex,SSSVertex>
describeHerwigNMSSMHHHVertex("Herwig::NMSSMHHHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMHHHVertex::Init() {

  static ClassDocumentation<NMSSMHHHVertex> documentation
    ("This is the triple Higgs coupling in the NMSSM.");

  static Switch<NMSSMHHHVertex,bool> interfaceIncludeRadiativeCorrections
    ("IncludeRadiativeCorrections",
     "Include radiative corrections in the vertex",
     &NMSSMHHHVertex::_includeRadiative, false, false, false);
  static SwitchOption interfaceIncludeRadiativeCorrectionsYes
    (interfaceIncludeRadiativeCorrections,
     "Yes",
     "Include the radiative terms",
     true);
  static SwitchOption interfaceIncludeRadiativeCorrectionsNo
    (interfaceIncludeRadiativeCorrections,
     "No",
     "Don't include them",
     false);

}

//calulate couplings
void NMSSMHHHVertex::setCoupling(Energy2 q2,tcPDPtr p1,tcPDPtr p2, 
				 tcPDPtr p3) {
  using Constants::pi;
  long higgs[3] = {p1->id(), p2->id(), p3->id()};
  unsigned int ns(0), np(0), nc(0);
  for( int i = 0; i < 3; ++i ) {
    if( higgs[i] == 25 || higgs[i] == 35 || higgs[i] == 45 )
      ++ns;
    else if( higgs[i] == 36 || higgs[i] == 46 )
      ++np;
    else if( abs(higgs[i]) ==  37 )
      ++nc;
  }
  //check three Higgs in vertex
  assert( ns + np + nc == 3 );
  if( q2 != _q2last ) {
    _q2last = q2;
    _glast = weakCoupling(q2);
    _mb = _theSM->mass(q2,getParticleData(5));
    _mt = _theSM->mass(q2,getParticleData(6));
  }
  //define VEV's
  double rt = sqrt(0.5);
  Energy _mtpole = getParticleData(6)->mass();
  Energy2 Qstsb = _MQ3*_MU2;
  double radlog(0.);
  if(_includeRadiative) {
    radlog = Qstsb/sqr(_mtpole);
    assert(radlog!=0.);
    radlog = log(radlog);
  }
  complex<Energy> coupling;
  //CP even Higgs
  if( ns == 3 ) {
    unsigned int a = (higgs[0] - 25)/10;
    unsigned int b = (higgs[1] - 25)/10;
    unsigned int c = (higgs[2] - 25)/10;
    coupling = 
      sqr(_lambda)*rt*(_vu*(usMix(a,b,c,1,0,0) + usMix(a,b,c,1,2,2))/_glast +
		       _vd*(usMix(a,b,c,0,1,1) + usMix(a,b,c,0,2,2))/_glast +
		       _s *(usMix(a,b,c,2,1,1) + usMix(a,b,c,2,0,0)))
      - _lambda*_kappa*rt*(_vu*usMix(a,b,c,0,2,2)/_glast +
			   _vd*usMix(a,b,c,2,1,2)/_glast + 2.*_s*usMix(a,b,c,1,0,2))
      + sqr(_kappa)/rt*_s*usMix(a,b,c,2,2,2)
      - _lambda*_theAl*rt*usMix(a,b,c,1,0,2)
      + _kappa*_theAk*rt/3.*usMix(a,b,c,2,2,2)
      + sqr(_glast)*0.25*rt/sqr(_cw)*(_vu*(usMix(a,b,c,1,1,1) -
					   usMix(a,b,c,1,0,0))/_glast -
				      _vd*(usMix(a,b,c,0,1,1) -
					   usMix(a,b,c,0,0,0))/_glast);
    // additional radiative terms
    if(_includeRadiative) {
      complex <Energy> radtop = usMix(a,b,c,1,1,1)*3.0*sqrt(2.0)*radlog
	*sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot= usMix(a,b,c,0,0,0)*3.0*sqrt(2.0)*radlog
	*sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      coupling += radbot + radtop;
    }					
  }
  //CP even, CP odd Vertex
  else if(ns == 1 && np == 2) {
    unsigned int a(0), b(0), c(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 ) {
      a = (higgs[0] - 25)/10;
      b = (higgs[1] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 ) {
      a = (higgs[1] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else {
      a = (higgs[2] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[1] - 36)/10;
    }
    coupling =	
      sqr(_lambda)*rt*(_vu*(upMix(a,b,c,1,0,0) + upMix(a,b,c,1,2,2))/_glast +
		       _vd*(upMix(a,b,c,0,1,1) + upMix(a,b,c,0,2,2))/_glast +
		       _s *(upMix(a,b,c,2,1,1) + upMix(a,b,c,2,0,0)))
      + _lambda*_kappa*rt*(_vu*(upMix(a,b,c,0,2,2) 
				- 2.*upMix(a,b,c,2,0,2))/_glast +
			   _vd*(upMix(a,b,c,1,2,2)
				- 2.*upMix(a,b,c,2,1,2))/_glast 
			   + 2.*_s*(upMix(a,b,c,2,1,0)
				    - upMix(a,b,c,1,0,2) - upMix(a,b,c,0,1,2)))
      + sqr(_kappa)/rt*_s*upMix(a,b,c,2,2,2)
      +_lambda*_theAl*rt*(upMix(a,b,c,1,0,2)
			  + upMix(a,b,c,0,1,2) + upMix(a,b,c,2,1,0))
      - _kappa*_theAk*rt*upMix(a,b,c,2,2,2)
      + sqr(_glast)*0.25*rt/sqr(_cw)*(_vu*(upMix(a,b,c,1,1,1) -
					   upMix(a,b,c,1,0,0))/_glast - 
				      _vd*(upMix(a,b,c,0,1,1) -
					   upMix(a,b,c,0,0,0))/_glast);
    if(_includeRadiative) {
      complex <Energy> radtop = upMix(a,b,c,1,1,1)*3.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot= upMix(a,b,c,0,0,0)*3.0*sqrt(2.0)*radlog*
	sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      coupling += radbot + radtop;
    }
  }
  //Charged Higgs
  else {
    unsigned int a(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 )
      a = (higgs[0] - 25)/10;
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 )
      a = (higgs[1] - 25)/10;
    else
      a = (higgs[2] - 25)/10;
    coupling = 	
      sqr(_lambda)*rt*2.*(_s*((*_mixS)(a,2)*sqr(_cb) + (*_mixS)(a,2)*sqr(_sb))
			  - (_vu*(*_mixS)(a,0)/_glast + 
			     _vd*(*_mixS)(a,1)/_glast)*_sb*_cb)
      +_lambda*_sb*_cb*2.*(*_mixS)(a,2)*(_kappa*_s/rt + rt*_theAl)
      + sqr(_glast)*0.5*rt*_sw2/sqr(_cw)*((_vu*(*_mixS)(a,1)/_glast - 
					   _vd*(*_mixS)(a,0)/_glast)*sqr(_cb) + 
					  (_vd*(*_mixS)(a,0)/_glast -
					   _vu*(*_mixS)(a,1)/_glast)*sqr(_sb))
      + sqr(_glast)*0.5*rt*(_vu*((*_mixS)(a,1)*sqr(_cb) + 
				 (*_mixS)(a,1)*sqr(_sb) + 
				 2.*(*_mixS)(a,0)*_cb*_sb)/_glast 
			    + _vd*((*_mixS)(a,0)*sqr(_cb) + 
				   (*_mixS)(a,0)*sqr(_sb) 
				   + 2.*(*_mixS)(a,1)*_sb*_cb)/_glast);
    if(_includeRadiative) {
      complex <Energy> radtop =(*_mixS)(a,1)*sqr(_sb)*6.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot=(*_mixS)(a,0)*sqr(_cb)*6.0*sqrt(2.0)*radlog*
	sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      complex <Energy> temp2 = _vu*((*_mixS)(a,1)*sqr(_cb) + 
				    (*_mixS)(a,0)*_sb*_cb)/_glast+ 
	_vd*((*_mixS)(a,1)*_sb*_cb + 
	     (*_mixS)(a,0)*sqr(_sb))/_glast;				
      
      complex <Energy> radtopbot= temp2*6.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mb)*sqr(_glast)*sqr(_glast)
	/(16.0*sqr(pi)*sqr(_vu)*sqr(_vd));
      coupling += radbot + radtop + radtopbot;					
    }
  }
  norm(-coupling * UnitRemoval::InvE);
}
