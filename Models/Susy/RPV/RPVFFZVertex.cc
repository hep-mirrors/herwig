// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFZVertex class.
//

#include "RPVFFZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "RPVhelper.h"

using namespace Herwig;

RPVFFZVertex::RPVFFZVertex()  : _sw(0.), _cw(0.), _id1last(0), 
				_id2last(0), _q2last(), _couplast(0.),
				_leftlast(0.), _rightlast(0.),
				_gl(17,0.0), _gr(17,0.0), _gblast(0),
				_interactions(0) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr RPVFFZVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFZVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << _theU << _theV 
     << _gl << _gr << _interactions;
}

void RPVFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> _theU >> _theV 
     >> _gl >> _gr >> _interactions;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVFFZVertex,Helicity::FFVVertex>
describeHerwigRPVFFZVertex("Herwig::RPVFFZVertex", "HwSusy.so HwRPV.so");

void RPVFFZVertex::Init() {

  static ClassDocumentation<RPVFFZVertex> documentation
    ("The RPVFFZVertex class implements trhe coupling of the Z to all"
     " fermion-antifermion pairs in models with bilinear RPV.");

  static Switch<RPVFFZVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVFFZVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsSM
    (interfaceInteractions,
     "SM",
     "Only include what would have been the interactions with the SM"
     " fermions in the absence of mixing",
     1);
  static SwitchOption interfaceInteractionsNeutralino
    (interfaceInteractions,
     "Neutralino",
     "Only include what would have been the interactions with the "
     "neutralinos in the absence of mixing",
     2);
  static SwitchOption interfaceInteractionsChargino
    (interfaceInteractions,
     "Chargino",
     "Only include what would have been the interactions with the "
     "charginos in the absence of mixing",
     3);
}

void RPVFFZVertex::doinit() {
  // extract the mixing matrices
  tSusyBasePtr model = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!model) throw InitException() << "RPVFFZVertex::doinit() - "
				   << "The model pointer is null."
				   << Exception::abortnow;
  _theN = model->neutralinoMix();
  _theU = model->charginoUMix();
  _theV = model->charginoVMix();
  if( !_theN || !_theU || !_theV )
    throw InitException() << "RPVFFZVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _theU << "  V: " << _theV << "  N: " << _theN
			  << Exception::abortnow;
  // Standard Model fermions
  if(_interactions==0||_interactions==1) {
    // PDG codes for the particles
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix, ix, 23);
    }
    // the leptons
    for(int ix=11;ix<17;ix+=2) {
      addToList(-ix, ix, 23);
    }
    for(int ix=12;ix<17;ix+=2) {
      if(_theN->size().first==7) {
	long inu = (ix-12)/2+17;
	addToList( inu, inu, 23);
      }
      else
	addToList(-ix, ix, 23);
    }
  }
  // neutralinos
  if(_interactions==0||_interactions==2) {
    vector<long> neu(4);
    neu[0] =  1000022; neu[1] = 1000023;
    neu[2] =  1000025; neu[3] = 1000035;
    if(_theN->size().first==7) {
      if(model->majoranaNeutrinos()) {
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
    for(unsigned int i = 0; i < neu.size(); ++i) {
      for(unsigned int j = 0; j < neu.size(); ++j) {
	if(!(i>3&&i==j)) addToList(neu[i], neu[j], 23);
      }
    }
  }
  // charginos
  if(_interactions==0||_interactions==3) {
    addToList(-1000024, 1000024, 22);
    addToList(-1000037, 1000037, 22);
    vector<long> cha(2);
    cha[0] = 1000024; cha[1] = 1000037;
    if(_theV->size().first==5) {
      cha.push_back(-11);
      cha.push_back(-13);
      cha.push_back(-15);
    }
    for(unsigned int i = 0; i < cha.size(); ++i) {
      for(unsigned int j = 0; j < cha.size(); ++j) {
	if(!(i>1&&i==j)) addToList(-cha[i], cha[j], 23);
      }
    }
  }
  Helicity::FFVVertex::doinit();
  // weak mixing
  double sw2 = sin2ThetaW();
  _cw  = sqrt(1. - sw2);
  _sw  = sqrt(   sw2  );
  // Standard Model couplings
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = -0.25*(model->vd()  + model->ad() );
    _gl[2*ix ]   = -0.25*(model->vu()  + model->au() );
    _gl[2*ix+9 ] = -0.25*(model->ve()  + model->ae() );
    _gl[2*ix+10] = -0.25*(model->vnu() + model->anu());
    _gr[2*ix-1]  = -0.25*(model->vd()  - model->ad() );
    _gr[2*ix ]   = -0.25*(model->vu()  - model->au() );
    _gr[2*ix+9 ] = -0.25*(model->ve()  - model->ae() );
    _gr[2*ix+10] = -0.25*(model->vnu() - model->anu());
  }
}

void RPVFFZVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  long iferm1(part1->id()), iferm2(part2->id()), boson(part3->id());
  long iferm = abs(iferm1);
  // chargino coupling to the photon
  if(part3->id()==ParticleID::gamma) {
    assert(iferm == abs(iferm2));
    _gblast = boson;
    _id1last = iferm1;
    _id2last = iferm2;
    _leftlast  = -1.;
    _rightlast = -1.;
    if(iferm1>0) {
      Complex temp = _leftlast;
      _leftlast  = -_rightlast;
      _rightlast = -temp;
    }
  }
  // coupling to the Z
  else {
    assert(part3->id()==ParticleID::Z0);
    // quarks
    if(iferm<=6) {
      _leftlast  = _gl[iferm]/(_sw*_cw);
      _rightlast = _gr[iferm]/(_sw*_cw);
    }
    // charged leptons and charginos
    else if(part1->iCharge()!=0) {
      if(boson != _gblast || iferm1 != _id1last || iferm2 != _id2last) {
	_gblast = boson;
	_id1last = iferm1;
	_id2last = iferm2;
	unsigned int ic1(0);
	if(_theV->size().first==2&&iferm<=16) {
	  _leftlast  = -_gr[iferm];
	  _rightlast = -_gl[iferm];
	}
	else {
	  ic1 = RPV_helper::charginoIndex(iferm1);
	  unsigned int ic2 = RPV_helper::charginoIndex(iferm2);
	  _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	    0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
	  _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	    0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
	  if(abs(iferm1) == abs(iferm2)) {
	    _leftlast  += sqr(_sw);
	    _rightlast += sqr(_sw);
	  }
	  if(_theV->size().first==5) {
	    for(unsigned int ix=0;ix<3;++ix) {
	      _rightlast += -0.5*(*_theU)(ic1, 2+ix)*conj((*_theU)(ic2, 2+ix));
	    }
	  }
	}
	if((ic1<2&&iferm1>0)||(ic1>=2&&iferm1<0)) {
	  Complex temp = _leftlast;
	  _leftlast  = -_rightlast;
	  _rightlast = -temp;
	}
	Complex temp = _leftlast;
	_leftlast  = -_rightlast;
	_rightlast = -temp;
	_leftlast  /= _sw*_cw;
	_rightlast /= _sw*_cw;
      }
    }
    // neutrinos and neutralinos
    else {
      // case where only 4x4 matrix and neutrino
      if(_theN->size().first==4&&iferm<=16) {
	assert(iferm==12||iferm==14||iferm==16);
	_leftlast  = _gl[iferm]/(_sw*_cw);
	_rightlast = _gr[iferm]/(_sw*_cw);
      }
      // neutralino
      else {
	long ic1 = part2->id();
	long ic2 = part1->id();
	assert(ic1 == ParticleID::SUSY_chi_10 || ic1 == ParticleID::SUSY_chi_20 ||
	       ic1 == ParticleID::SUSY_chi_30 || ic1 == ParticleID::SUSY_chi_40 ||
	       abs(ic1) == 12 || abs(ic1) == 14 || abs(ic1) == 16 ||
	       abs(ic1) == 17 || abs(ic1) == 18 || abs(ic1) == 19 );
	assert(ic2 == ParticleID::SUSY_chi_10 || ic2 == ParticleID::SUSY_chi_20 ||
	       ic2 == ParticleID::SUSY_chi_30 || ic2 == ParticleID::SUSY_chi_40 ||
	       abs(ic2) == 12 || abs(ic2) == 14 || abs(ic2) == 16 ||
	       abs(ic2) == 17 || abs(ic2) == 18 || abs(ic2) == 19 );
	if(ic1 != _id1last || ic2 != _id2last) {
	  _id1last = ic1;
	  _id2last = ic2;
	  unsigned int neu1 = RPV_helper::neutralinoIndex(ic1);
	  unsigned int neu2 = RPV_helper::neutralinoIndex(ic2);
	  _leftlast = 0.5*( (*_theN)(neu1, 3)*conj((*_theN)(neu2, 3)) -
			    (*_theN)(neu1, 2)*conj((*_theN)(neu2, 2)) );
	  if(_theN->size().first>4) {
	    for(unsigned int k=0;k<3;++k)
	      _leftlast -= 0.5*(*_theN)(neu1, 4+k)*conj((*_theN)(neu2, 4+k));
	  }
	  _rightlast = -conj(_leftlast);
	  _leftlast  /= _sw*_cw;
	  _rightlast /= _sw*_cw;
	}
      }
    }
  }
  norm ( _couplast);
  left ( _leftlast);
  right(_rightlast);
}
