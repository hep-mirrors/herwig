// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFSVertex class.
//

#include "RPVFFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {

  unsigned int neutralinoIndex(long id) {
    if(id> 1000000)
      return id<1000025 ? id-1000022 : (id-1000005)/10;
    else if(abs(id)<=16) 
      return (abs(id)-4)/2;
    else
      return id-13;
  }
  
  unsigned int charginoIndex(long id) {
    return abs(id)>1000000 ? (abs(id)-1000024)/13 : (abs(id)-7)/2;
  }

}

RPVFFSVertex::RPVFFSVertex() : interactions_(0), mw_(ZERO),
			       _q2last(ZERO), _couplast(0.),
			       _leftlast(0.),_rightlast(0.),
			       _id1last(0), _id2last(0), _id3last(0),
			       yukawa_(true),
			       _sw(0.), _cw(0.),_sb(0.), _cb(0.),
			       vd_(ZERO), vu_(ZERO),
			       _massLast(make_pair(ZERO,ZERO)) {
  orderInGem(1);
  orderInGs(0);
}  

IBPtr RPVFFSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFSVertex::doinit() {
  // cast the model to the RPV one
  model_ = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model_ ) throw InitException() << "RPVFFSVertex::doinit() - "
				      << "The pointer to the MSSM object is null!"
				      << Exception::abortnow;
  // get the various mixing matrices
  _stop = model_->stopMix();
  _sbot = model_->sbottomMix();
  _stau = model_->stauMix();
  _nmix = model_->neutralinoMix();
  _umix = model_->charginoUMix();
  _vmix = model_->charginoVMix();
  _mixH = model_->CPevenHiggsMix();
  _mixP = model_->CPoddHiggsMix();
  _mixC = model_->ChargedHiggsMix();
  if(!_stop || !_sbot ) throw InitException() << "RPVFFSVertex::doinit() - "
					      << "A squark mixing matrix pointer is null."
					      << " stop: " << _stop << " sbottom: "
					      << _sbot << Exception::abortnow;
  if(!_nmix) throw InitException() << "RPVFFSVertex::doinit() - "
				   << "A gaugino mixing matrix pointer is null."
				   << " N: " << _nmix << " U: " << _umix 
				   << " V: " << _vmix << Exception::abortnow;
  if(! ( _stau || (!_stau&& _mixC->size().first>1))) 
    throw InitException() << "RPVFFSVertex::doinit() - "
			  << "Must have either the stau mixing matrix or it"
			  << " must be part of the charged Higgs boson mixing."
			  << Exception::abortnow;
  // various interactions
  // scalar Higgs bosons
  vector<long> h0(2);
  h0[0] = 25; h0[1] = 35; 
  if(_mixH&&_mixH->size().first>2) {
    h0.push_back(1000012); h0.push_back(1000014); h0.push_back(1000016);
  }
  // pseudoscalar Higgs bosons
  vector<long> a0(1,36);
  if(_mixP&&_mixP->size().first>1) {
    a0.push_back(1000017); a0.push_back(1000018); a0.push_back(1000019);
  }
  // charged Higgs bosons
  vector<long> hp(1,37);
  if(_mixC->size().first>1) {
    hp.push_back(-1000011); hp.push_back(-1000013); hp.push_back(-1000015);
    hp.push_back(-2000011); hp.push_back(-2000013); hp.push_back(-2000015);
  }
  // neutralinos
  vector<long> neu(4);
  neu[0] = 1000022; neu[1] = 1000023;
  neu[2] = 1000025; neu[3] = 1000035;
  if(_nmix->size().first>4) {
    if(model_->majoranaNeutrinos()) {
      neu.push_back(17);
      neu.push_back(18);
      neu.push_back(19);
    }
    else {
      neu.push_back(12);
      neu.push_back(14);
      neu.push_back(16);
    }
  }
  // charginos
  vector<long> chg(2);
  chg[0] = 1000024; chg[1] = 1000037;
  if(_umix->size().first>2) {
    chg.push_back(-11); chg.push_back(-13); chg.push_back(-15);
  }
  // FFH
  if(interactions_==0||interactions_==1) {
    // quarks neutral scalar
    for ( unsigned int h = 0; h < h0.size(); ++h ) {
      for(long ix=1;ix<7;++ix) addToList(-ix,ix,h0[h]);
    } 
    // quarks neutral pseudoscalar
    for ( unsigned int h = 0; h < a0.size(); ++h ) {
      for(long ix=1;ix<7;++ix) addToList(-ix,ix,a0[h]);
    }
    // quarks charged higgs
    for(unsigned int h=0;h<hp.size();++h) {
      for(long ix=1;ix<6;ix+=2) {
	//outgoing H+
	addToList(-ix-1,  ix, hp[h]);
	//outgoing H-
	addToList(-ix  ,ix+1,-hp[h]);
      }
    }
    // charged leptons neutral scalar
    for ( unsigned int h = 0; h < h0.size(); ++h ) {
      for(long ix=11;ix<16;ix+=2) addToList(-ix,ix,h0[h]);
    } 
    // charged leptons neutral pseudoscalar
    for ( unsigned int h = 0; h < a0.size(); ++h ) {
      for(long ix=11;ix<16;ix+=2) addToList(-ix,ix,a0[h]);
    }
    // charged higgs to leptons, no mixing
    if(_nmix->size().first<=4) {
      for(unsigned int h=0;h<hp.size();++h) {
	for(long ix=11;ix<16;ix+=2) {
	  //outgoing H+
	  addToList(-ix-1,  ix, hp[h]);
	  //outgoing H-
	  addToList(-ix  ,ix+1,-hp[h]);
	}
      }
    }
  }
  // GOGOH
  if(interactions_==0||interactions_==2) {
    // neutral scalar Higgs neutralino
    for(unsigned int i=0;i<h0.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = j; k < neu.size(); ++k) {
	  addToList(neu[j], neu[k], h0[i]);
	}
      }
    }
    // neutral pseudoscalar Higgs neutralino
    for(unsigned int i=0;i<a0.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = j; k < neu.size(); ++k) {
	  addToList(neu[j], neu[k], a0[i]);
	}
      }
    }
    // neutral scalar Higgs chargino
    for(unsigned int i=0;i<h0.size();++i) {
      for(unsigned int j = 0; j < chg.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
	  if(j==k&&abs(chg[j])<16) continue;
	  addToList(-chg[j], chg[k], h0[i]);
	}
      }
    }
    // neutral scalar Higgs chargino
    for(unsigned int i=0;i<a0.size();++i) {
      for(unsigned int j = 0; j < chg.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
	  if(j==k&&abs(chg[j])<16) continue;
	  addToList(-chg[j], chg[k], a0[i]);
	}
      }
    }
    // charged higgs
    for(unsigned int i=0;i<hp.size();++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	for(unsigned int k = 0; k < chg.size(); ++k) {
 	  addToList(-chg[k], neu[j], hp[i]);
 	  addToList( chg[k], neu[j],-hp[i]);
	}
      }
    }
  }
  // neutralino sfermion
  if(interactions_==0||interactions_==3) {
    // quarks
    for(unsigned int nl = 0; nl < neu.size(); ++nl) {
      for(long ix=1;ix<7;++ix){
	addToList( neu[nl],  ix, -(1000000+ix) );
	addToList( neu[nl],  ix, -(2000000+ix) );
	addToList( neu[nl], -ix,  (1000000+ix) );
	addToList( neu[nl], -ix,  (2000000+ix) );
      }
    }
    // neutrino
    if(_nmix->size().first<=4) {
      for(unsigned int nl = 0; nl < neu.size(); ++nl) {
	for(long ix=12;ix<17;ix+=2) {
	  addToList( neu[nl],  ix, -(1000000+ix) );
	  addToList( neu[nl], -ix,  (1000000+ix) );
	}
      }
    }
    // charged leptons no mixing
    if(!_mixC || _mixC->size().first==1) {
      for(unsigned int nl = 0; nl < neu.size(); ++nl) {
	for(long ix=11;ix<17;ix+=2) {
	  addToList( neu[nl],  ix, -(1000000+ix) );
	  addToList( neu[nl], -ix,  (1000000+ix) );
	  addToList( neu[nl],  ix, -(2000000+ix) );
	  addToList( neu[nl], -ix,  (2000000+ix) );
	}
      }
    }
  }
  // chargino sfermion
  if(interactions_==0||interactions_==4) {
    //quarks 
    for(unsigned int ic = 0; ic < chg.size(); ++ic) {
      for(long ix = 1; ix < 7; ++ix) {
	if( ix % 2 == 0 ) {
	  addToList(-chg[ic], ix,-( 999999+ix));
	  addToList(-chg[ic], ix,-(1999999+ix));
	  addToList( chg[ic],-ix,   999999+ix );
	  addToList( chg[ic],-ix,  1999999+ix );
	}
	else {
	  addToList(-chg[ic],-ix,  1000001+ix );
	  addToList(-chg[ic],-ix,  2000001+ix );
	  addToList( chg[ic], ix,-(1000001+ix));
	  addToList( chg[ic], ix,-(2000001+ix));
	}
      }
    }
    // sneutrinos
    if(!_mixH || _mixH->size().first<=2) {
      for(unsigned int ic = 0; ic < chg.size(); ++ic) {
	for(long ix = 11; ix < 17; ix+=2) {
	  addToList(-chg[ic],-ix,1000001+ix);
	  addToList(chg[ic],ix,-(1000001+ix));
	}
      }
    }
    // charged leptons
    if(!_mixC || _mixC->size().first==1) {
      for(unsigned int ic = 0; ic < chg.size(); ++ic) {
	for(long ix = 12; ix < 17; ix+=2) {
	  addToList(-chg[ic], ix,-( 999999+ix));
	  addToList(-chg[ic], ix,-(1999999+ix));
	  addToList( chg[ic],-ix, ( 999999+ix));
	  addToList( chg[ic],-ix, (1999999+ix));
	}
      }
    }
  }
  FFSVertex::doinit();
  // various couplings and parameters
  mw_ = getParticleData(ParticleID::Wplus)->mass();
  _sw = sqrt(sin2ThetaW());
  double tb = model_->tanBeta();
  _cw = sqrt(1. - sqr(_sw));
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  vector<Energy> vnu = model_->sneutrinoVEVs();
  double g = electroMagneticCoupling(sqr(mw_))/_sw;
  Energy v = 2.*mw_/g;
  vd_ = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tb)));
  vu_ = vd_*tb;
  // couplings of the neutral scalar Higgs to charginos
  Energy me   = model_->mass(sqr(mw_),getParticleData(ParticleID::eminus  ));
  Energy mmu  = model_->mass(sqr(mw_),getParticleData(ParticleID::muminus ));
  Energy mtau = model_->mass(sqr(mw_),getParticleData(ParticleID::tauminus));
  double h_E[3] = {sqrt(2.)*me/vd_/g,sqrt(2.)*mmu /vd_/g,sqrt(2.)*mtau/vd_/g};
  for(unsigned int ih=0;ih<_mixH->size().first;++ih ) {
    OCCHL_.push_back(vector<vector<Complex> >
		     (_umix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i=0;i<_umix->size().first;++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCCHL_[ih][i][j]    = -sqrt(0.5)*
	  ( (*_mixH)(ih,0)*(*_vmix)(i,0)*(*_umix)(j,1) +
	    (*_mixH)(ih,1)*(*_vmix)(i,1)*(*_umix)(j,0));
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCCHL_[ih][i][j] -= sqrt(0.5)*(+(*_mixH)(ih,k)*(*_vmix)(i,0)*(*_umix)(j,k)
					 +h_E[k-2]*(*_umix)(j,k)*(*_vmix)(i,k)*(*_mixH)(ih,0)
					 -h_E[k-2]*(*_umix)(j,1)*(*_vmix)(i,k)*(*_mixH)(ih,k));
	}
      }
    }
  }
  // couplings of the neutral scalar Higgs to neutralinos
  double tw = _sw/_cw;
  for(unsigned int ih=0;ih<_mixH->size().first;++ih) {
    ONNHL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_nmix->size().first,0.)));
    for(unsigned int i=0;i<_nmix->size().first;++i) {
      for(unsigned int j=0;j<_nmix->size().first;++j) {
	ONNHL_[ih][i][j] =0.;
	for(unsigned int in=2;in<_nmix->size().first;++in) {
	  double sign = in!=3 ? 1. : -1.;
	  ONNHL_[ih][i][j] += 0.5*sign*(*_mixH)(ih,in-2)*
	    ( + (tw*(*_nmix)(i,0) - (*_nmix)(i,1) )*(*_nmix)(j,in)
	      + (tw*(*_nmix)(j,0) - (*_nmix)(j,1) )*(*_nmix)(i,in) );
	} 
      }
    }
  }
  // couplings of the neutral pseudoscalar Higgs to neutralinos
  for(unsigned int ih=0;ih<_mixP->size().first;++ih) {
    ONNAL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_nmix->size().first,0.)));
    for(unsigned int i=0;i<_nmix->size().first;++i) {
      for(unsigned int j=0;j<_nmix->size().first;++j) {
	ONNAL_[ih][i][j] =0.;
	for(unsigned int in=2;in<_nmix->size().first;++in) {
	  double sign = in!=3 ? 1. : -1.;
	  ONNAL_[ih][i][j] += -0.5*sign*(*_mixP)(ih,in-2)*
	    ( + (tw*(*_nmix)(i,0) - (*_nmix)(i,1) )*(*_nmix)(j,in)
	      + (tw*(*_nmix)(j,0) - (*_nmix)(j,1) )*(*_nmix)(i,in) );
	}
      }
    }
  }
  // couplings of the neutral pseudoscalar higgs to charginos
  for(unsigned int ih=0;ih<_mixP->size().first;++ih) {
    OCCAL_.push_back(vector<vector<Complex> >
		     (_umix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i=0;i<_umix->size().first;++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCCAL_[ih][i][j] = 
	  (*_mixP)(ih,0)*(*_vmix)(i,0)*(*_umix)(j,1)+
	  (*_mixP)(ih,1)*(*_vmix)(i,1)*(*_umix)(j,0);
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCCAL_[ih][i][j] += (*_mixP)(ih,k)*(*_vmix)(i,0)*(*_umix)(j,k)
	    -h_E[k-2]*(*_umix)(j,k)*(*_vmix)(i,k)*(*_mixP)(ih,0)
	    +h_E[k-2]*(*_umix)(j,1)*(*_vmix)(i,k)*(*_mixP)(ih,k);
	}
	OCCAL_[ih][i][j] *= sqrt(0.5);
      }
    }
  }
  // couplings for the charged higgs
  for(unsigned int ih=0;ih<_mixC->size().first;++ih) {
    OCNSL_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_umix->size().first,0.)));
    OCNSR_.push_back(vector<vector<Complex> >
		     (_nmix->size().first,vector<Complex>(_umix->size().first,0.)));
    for(unsigned int i = 0; i < _nmix->size().first; ++i) {
      for(unsigned int j=0;j<_umix->size().first;++j) {
	OCNSL_[ih][i][j] = (*_mixC)(ih,1)*conj((*_nmix)(i, 3)*(*_vmix)(j,0)
					       +((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_vmix)(j,1)/sqrt(2));


	OCNSR_[ih][i][j] = (*_mixC)(ih,0)*    ((*_nmix)(i, 2)*(*_umix)(j,0)
					       -((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_umix)(j,1)/sqrt(2));
	for(unsigned int k=2;k<_umix->size().first;++k) {
	  OCNSL_[ih][i][j] += -h_E[k-2]*(*_mixC)(ih,0)*conj((*_nmix)(i,2+k)*(*_vmix)(j,k))
	    +(*_mixC)(ih,k  )*h_E[k-2]*conj((*_nmix)(i, 2)*(*_vmix)(j,k))
	    +(*_mixC)(ih,k+3)*tw*sqrt(2.)*conj((*_nmix)(i,0)*(*_vmix)(j,k));
	  OCNSR_[ih][i][j] += (*_mixC)(ih,k)*((*_nmix)(i,2+k)*(*_umix)(j,0)
					      -((*_nmix)(i,1) + (*_nmix)(i,0)*tw)*(*_umix)(j,k)/sqrt(2))
	    -(*_mixC)(ih,k+3)*h_E[k-2]*((*_nmix)(i,k+2)*(*_umix)(j,1)-(*_nmix)(i,2)*(*_umix)(j,k));
	}
      }
    }
  }
}

void RPVFFSVertex::persistentOutput(PersistentOStream & os) const {
  os << interactions_ << _stop << _sbot << _stau << _umix << _vmix
     << _nmix << _mixH << _mixP << _mixC << ounit(mw_,GeV) << yukawa_
     << model_  << _sw << _cw << _sb << _cb << ounit(vd_,GeV) << ounit(vu_,GeV)
     << OCCHL_ << ONNHL_ << ONNAL_ << OCCAL_ << OCNSL_ << OCNSR_;
}

void RPVFFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> interactions_ >> _stop >> _sbot >> _stau >> _umix >> _vmix
     >> _nmix >> _mixH >> _mixP >> _mixC >> iunit(mw_,GeV) >> yukawa_
     >> model_ >> _sw >> _cw  >> _sb >> _cb >> iunit(vd_,GeV) >> iunit(vu_,GeV)
     >> OCCHL_ >> ONNHL_ >> ONNAL_ >> OCCAL_ >> OCNSL_ >> OCNSR_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVFFSVertex,Helicity::FFSVertex>
describeHerwigRPVFFSVertex("Herwig::RPVFFSVertex", "HwSusy.so HwRPV.so");

void RPVFFSVertex::Init() {

  static ClassDocumentation<RPVFFSVertex> documentation
    ("The RPVFFSVertex class implements all the couplings of fermion-antiferion to scalars"
     " in R-Parity violating models, including sfermion fermion gaugino, SM ferimon antiferimon Higgs "
     "and gaugino-gaugino Higgs.");

  static Switch<RPVFFSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Whice interactions to include",
     &RPVFFSVertex::interactions_, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsHiggsSMFermions
    (interfaceInteractions,
     "HiggsSMFermions",
     "Interactions of Higgs with SM fermions",
     1);
  static SwitchOption interfaceInteractionsHiggsGaugino
    (interfaceInteractions,
     "HiggsGaugino",
     "Interactions of the Higgs with the gauginos",
     2);
  static SwitchOption interfaceInteractionsNeutralinoSfermion
    (interfaceInteractions,
     "NeutralinoSfermion",
     "Include the neutralino sfermion interactions",
     3);
  static SwitchOption interfaceInteractionsCharginoSfermions
    (interfaceInteractions,
     "CharginoSfermions",
     "Include the chargino sfermion interactions",
     4);

  static Switch<RPVFFSVertex,bool> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings in neutralino/chargino interactions",
     &RPVFFSVertex::yukawa_, true, false, false);
  static SwitchOption interfaceYukawaYes
    (interfaceYukawa,
     "Yes",
     "Include the terms",
     true);
  static SwitchOption interfaceYukawaNo
    (interfaceYukawa,
     "No",
     "Don't include them",
     false);
}

void RPVFFSVertex::neutralinoSfermionCoupling(Energy2 q2, tcPDPtr fermion,
					      tcPDPtr gaugino, tcPDPtr sfermion) {
  long ism(abs(fermion->id())),ig(abs(gaugino->id())),isc(sfermion->id());
  if( ig != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ig;
    _id2last = ism;
    _id3last = isc;
    // sfermion mass eigenstate
    unsigned int alpha(abs(isc)/1000000 - 1);
    // neutralino state
    unsigned int nl = neutralinoIndex(ig);
    assert(nl<=6);
    // common primed neutralino matrices
    Complex n2prime = (*_nmix)(nl,1)*_cw - (*_nmix)(nl,0)*_sw;
    //handle neutrinos first
    if( ism == 12 || ism == 14 || ism == 16 ) {
      _leftlast = Complex(0., 0.);
      _rightlast = -sqrt(0.5)*n2prime/_cw;
    }
    else {
      Complex n1prime = (*_nmix)(nl,0)*_cw + (*_nmix)(nl,1)*_sw;
      tcPDPtr smf = getParticleData(ism);
      double qf = smf->charge()/eplus;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double y = yukawa_ ? double(model_->mass(q2, smf)/2./mw_) : 0.;
      double lambda(0.);
      //neutralino mixing element
      Complex nlf(0.);
      if( ism % 2 == 0 ) {
	y /= _sb;
	lambda = -0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,3);
      }
      else { 
	y /= _cb;
	lambda = 0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,2);
      }
      Complex bracketr = _sw*qf*n1prime - n2prime*lambda/_cw;
      //heavy quarks/sleptons
      if( ism == 5 || ism == 6 || ism == 15 ) {
	Complex ma1(0.), ma2(0.);
	if( ism == 5 ) {
	  ma1 = (*_sbot)(alpha,0);
	  ma2 = (*_sbot)(alpha,1);
	} 
	else if( ism == 6 ) {
	  ma1 = (*_stop)(alpha,0);
	  ma2 = (*_stop)(alpha,1);
	} 
	else {
	  ma1 = (*_stau)(alpha,0);
	  ma2 = (*_stau)(alpha,1);
	}
	_leftlast = y*conj(nlf)*ma1 - ma2*bracketl;
	_rightlast = y*nlf*ma2 + ma1*bracketr;
      }
      else {
	if( alpha == 0 ) {
	  _leftlast = y*conj(nlf);
	  _rightlast = bracketr;
	} 
	else {
	  _leftlast = -bracketl;
	  _rightlast = y*nlf;
	}
      }
      _leftlast  *= -sqrt(2.);
      _rightlast *= -sqrt(2.);
    }
  }
  //determine the helicity order of the vertex
  if( fermion->id() < 0 ) {
    left (conj(_rightlast));
    right(conj( _leftlast));
  }
  else {
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::charginoSfermionCoupling(Energy2 q2, tcPDPtr fermion,
					    tcPDPtr gaugino, tcPDPtr sfermion) {
  long ism(abs(fermion->id())),ig(abs(gaugino->id())),isc(sfermion->id());
  if( ig != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ig;
    _id2last = ism;
    _id3last = isc;
    // sfermion mass eigenstate
    unsigned int alpha(abs(isc)/1000000 - 1);
    // get the type of chargino
    unsigned int ch = charginoIndex(ig);
    assert(ch<=4);
    // various mixing matrix elements
    Complex ul1 = (*_umix)(ch,0), ul2 = (*_umix)(ch,1);
    Complex vl1 = (*_vmix)(ch,0), vl2 = (*_vmix)(ch,1);
    // lepton/slepton
    if( ism >= 11 && ism <= 16 ) {
      long lept = ( ism % 2 == 0 ) ? ism - 1 : ism;
      double y = yukawa_ ? 
	double(model_->mass(q2, getParticleData(lept))/mw_/sqrt(2)/_cb) : 0.;
      if( ism == 12 || ism == 14 ) {
	_leftlast = 0.;
	_rightlast = alpha == 0 ? - ul1 : y*ul2;
      }
      else if( ism == 16 ) {
	_leftlast = 0.;
	_rightlast = -ul1*(*_stau)(alpha, 0) + y*(*_stau)(alpha,1)*ul2;
      }
      else if( ism == 11 || ism == 13 || ism == 15 ) {
	_leftlast = y*conj(ul2);
	_rightlast = -vl1;
      }
    }
    // squark/quark
    else {
      double yd(0.), yu(0.);
      if(yukawa_) {
	if( ism % 2 == 0) {
	  yu = model_->mass(q2, getParticleData(ism))/mw_/sqrt(2)/_sb;
	  yd = model_->mass(q2, getParticleData(ism - 1))/mw_/sqrt(2)/_cb;
	}
	else {
	  yu = model_->mass(q2, getParticleData(ism + 1))/mw_/sqrt(2)/_sb;
	  yd = model_->mass(q2, getParticleData(ism))/mw_/sqrt(2)/_cb;
	}
      }
      //heavy quarks
      if( ism == 5 ) {
	_leftlast =  yd*conj(ul2)*(*_stop)(alpha,0);
	_rightlast = -vl1*(*_stop)(alpha, 0) + yu*vl2*(*_stop)(alpha,1);
      }
      else if( ism == 6 ) {
	_leftlast =  yu*conj(vl2)*(*_sbot)(alpha,0);
	_rightlast = -ul1*(*_sbot)(alpha, 0) + yd*ul2*(*_sbot)(alpha,1);
      }
      else {
	if( alpha == 0 ) {
	  _leftlast  = (ism % 2 == 0) ? yu*conj(vl2) : yd*conj(ul2);
	  _rightlast = (ism % 2 == 0) ? -ul1 : -vl1;
	}
	else {
	  _leftlast = 0.;
	  _rightlast = (ism % 2 == 0) ? yd*ul2 : yu*vl2;
	}
      }
    }
  }
  //determine the helicity order of the vertex
  if( fermion->id() < 0 ) {
    left (conj(_rightlast));
    right(conj( _leftlast));
  }
  else {
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::higgsFermionCoupling(Energy2 q2, tcPDPtr f1,
					tcPDPtr f2, tcPDPtr higgs) {
  long f1ID(f1->id()), f2ID(f2->id()), isc(higgs->id());
  // running fermion masses
  if( q2 != _q2last || _id1last  != f1ID) {
    _massLast.first  = model_->mass(q2,f1);
    _id1last  = f1ID;
  }
  if( q2 != _q2last || _id2last != f2ID) {
    _massLast.second = model_->mass(q2,f2);
    _id2last = f2ID;
  }
  if( q2 != _q2last) _id3last = isc;
  Complex fact(0.);
  // scalar neutral Higgs
  if(isc == ParticleID::h0         || isc  == ParticleID::H0         ||
     isc == ParticleID::SUSY_nu_eL || isc == ParticleID::SUSY_nu_muL ||
     isc == ParticleID::SUSY_nu_tauL ) {
    unsigned int ih = isc < 1000000 ? (isc-25)/10 : (isc-1000008)/2;
    unsigned int id = abs(f1ID);
    fact = -_massLast.first*
      ((id%2==0) ? (*_mixH)(ih,1)/vu_ : (*_mixH)(ih,0)/vd_);
    left (1.);
    right(1.);
  }
  // pseudoscalar neutral Higgs
  else if(isc == ParticleID::A0 || isc == 1000017 || isc == 1000018 ||
	  isc == 1000019 ) {
    unsigned int ih = isc < 1000000 ? 0 : (isc-1000016);
    unsigned int id = abs(f1ID);
    if(_mixP) {
      fact = -Complex(0., 1.)*_massLast.first*
	( (id%2==0) ?  (*_mixP)(ih,1)/vu_ : (*_mixP)(ih,0)/vd_);
    }
    else {
      fact = -Complex(0., 1.)*_massLast.first*
	( (id%2==0) ?  _cb/vu_ : _sb/vd_);
    }
    left ( 1.);
    right(-1.);
  }
  // charged higgs
  else {
    if(!_mixC) {
      if( abs(f1ID) % 2 == 0 ) {
	_leftlast  =  _massLast.first /vu_*_cb;
	_rightlast =  _massLast.second/vd_*_sb;
      }
      else {
	_leftlast  =  _massLast.second/vu_*_cb;
	_rightlast =  _massLast.first /vd_*_sb;
      }
    }
    else {
      unsigned int ih;
      if(abs(isc)==ParticleID::Hplus) {
	ih = 0;
      }
      else {
	isc *= -1;
	ih = abs(isc)<2000000 ? (abs(isc)-1000009)/2 : (abs(isc)-2000003)/2;
      }
      if( abs(f1ID) % 2 == 0 ) {
	_leftlast  =  _massLast.first /vu_*(*_mixC)(ih,1);
	_rightlast =  _massLast.second/vd_*(*_mixC)(ih,0);
      }
      else {
	_leftlast  =  _massLast.second/vu_*(*_mixC)(ih,1);
	_rightlast =  _massLast.first /vd_*(*_mixC)(ih,0);
      }
    }
    if( isc > 0 ) swap(_leftlast,_rightlast);
    fact = sqrt(2.);
    left ( _leftlast);
    right(_rightlast);
  }
  norm(fact);
}

void RPVFFSVertex::higgsGauginoCoupling(Energy2, tcPDPtr f1,
					tcPDPtr f2, tcPDPtr higgs) {
  long f1ID(f1->id()), f2ID(f2->id()), isc(higgs->id());
  if( isc == _id3last && f1ID == _id1last && f2ID == _id2last ) {
    left ( _leftlast);
    right(_rightlast);
  }
  else {
    _id1last = f1ID;
    _id2last = f2ID;
    _id3last = isc;
    // scalar neutral Higgs
    if(isc == ParticleID::h0         || isc  == ParticleID::H0         ||
       isc == ParticleID::SUSY_nu_eL || isc == ParticleID::SUSY_nu_muL ||
       isc == ParticleID::SUSY_nu_tauL ) {
      unsigned int ih = isc < 1000000 ? (isc-25)/10 : (isc-1000008)/2;
      // charginos
      if(f1->charged()) {
	unsigned int ei = charginoIndex(f1ID);
	unsigned int ej = charginoIndex(f2ID);
	if     (ei< 2&&f1ID>0) swap(ei,ej);
	else if(ei>=2&&f1ID<0) swap(ei,ej);
	_rightlast  = conj(OCCHL_[ih][ej][ei]);
	_leftlast   =      OCCHL_[ih][ei][ej] ;
      }
      // neutralinos
      else {
	unsigned int ei = neutralinoIndex(f1ID);
	unsigned int ej = neutralinoIndex(f2ID);
	_leftlast  = conj(ONNHL_[ih][ej][ei]);
	_rightlast =      ONNHL_[ih][ei][ej] ;
      }
    }
    // pseudoscalar neutral Higgs
    else if(isc == ParticleID::A0 || isc == 1000017 || isc == 1000018 ||
	    isc == 1000019 ) {
      unsigned int ih = isc < 1000000 ? 0 : (isc-1000016);
      // charginos
      if(f1->charged()) {
	unsigned int ei = charginoIndex(f1ID);
	unsigned int ej = charginoIndex(f2ID);
	if     (ei< 2&&f1ID>0) swap(ei,ej);
	else if(ei>=2&&f1ID<0) swap(ei,ej);
	_rightlast = -Complex(0.,1.)*conj(OCCAL_[ih][ej][ei]);
	_leftlast  =  Complex(0.,1.)*     OCCAL_[ih][ei][ej] ;
      }
      // neutralinos
      else {
	unsigned int ei = neutralinoIndex(f1ID);
	unsigned int ej = neutralinoIndex(f2ID);
	_leftlast  =  Complex(0.,1.)*conj(ONNAL_[ih][ej][ei]);
	_rightlast = -Complex(0.,1.)*     ONNAL_[ih][ei][ej] ;
      }
    }
    // charged higgs
    else {
      unsigned int ih = abs(isc) < 1000000 ? 0 : 
	(abs(isc) < 2000000 ? (abs(isc)-1000009)/2 : (abs(isc)-2000003)/2);
      long chg(f2ID), neu(f1ID);
      if(f1->charged()) swap(chg, neu);
      unsigned int ei = neutralinoIndex(neu);
      unsigned int ej = charginoIndex(chg);
      _leftlast  = -OCNSL_[ih][ei][ej];
      _rightlast = -OCNSR_[ih][ei][ej];
      bool chargedSwap = abs(isc)<1000000 ? isc<0 : isc>0;
      if( chargedSwap ) {
	Complex tmp = _leftlast;
	_leftlast  = conj(_rightlast);
	_rightlast = conj(tmp);
      }
    }
    left ( _leftlast);
    right(_rightlast);
  }
  norm(_couplast);
}

void RPVFFSVertex::setCoupling(Energy2 q2, tcPDPtr part1, 
			       tcPDPtr part2,tcPDPtr part3) {
  // overall normalisation
  if(q2!=_q2last || _couplast==0.) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  long f1ID(part1->id()), f2ID(part2->id()), isc(abs(part3->id()));
  // squark quark
  if(part3->coloured()) {
    tcPDPtr smfermion = part1, gaugino = part2;
    if(gaugino->coloured()) swap(smfermion,gaugino);
    if(gaugino->charged())
      charginoSfermionCoupling(q2,smfermion,gaugino,part3);
    else
      neutralinoSfermionCoupling(q2,smfermion,gaugino,part3);
  }
  // slepton/lepton without mixing
  else if((( isc >= 1000011 && isc <= 1000016) ||
	   ( isc >= 2000011 && isc <= 2000016)) && 
	  (!_mixC || _mixC->size().first<=1 || 
	   !_mixP || _mixP->size().first<=1 )) {
    tcPDPtr smfermion = part1, gaugino = part2;
    if(abs(gaugino->id())<1000000) swap(smfermion,gaugino);
    if(gaugino->charged())
      charginoSfermionCoupling(q2,smfermion,gaugino,part3);
    else
      neutralinoSfermionCoupling(q2,smfermion,gaugino,part3);
  }
  // SM quarks and Higgs
  else if((abs(f1ID) <=  6 && abs(f2ID) <=  6) ||
	  ((abs(f1ID) >= 11 && abs(f1ID) <= 16) &&
	   (abs(f2ID) >= 11 && abs(f2ID) <= 16) && 
	   _umix->size().first==2) ) {
    higgsFermionCoupling(q2,part1,part2,part3);
  }
  // gauginos and the Higgs (general case for sleptons)
  else {
    higgsGauginoCoupling(q2,part1,part2,part3);
  }
}
