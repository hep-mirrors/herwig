// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVFFWVertex class.
//

#include "RPVFFWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RPVhelper.h"

using namespace Herwig;

RPVFFWVertex::RPVFFWVertex() : _diagonal(false), _ckm(3,vector<Complex>(3,0.0)),
			       _sw(0.), _couplast(0.), _q2last(ZERO), 
			       _id1last(0), _id2last(0), _leftlast(0.),
			       _rightlast(0.), _interactions(0) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

IBPtr RPVFFWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVFFWVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVFFWVertex::doinit() {
  // SUSY mixing matrices
  tSusyBasePtr model = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "RPVFFWVertex::doinit() - The model pointer is null!"
			  << Exception::abortnow;
  _theN = model->neutralinoMix();
  _theU = model->charginoUMix();
  _theV = model->charginoVMix();
  if(!_theN || !_theU || ! _theV)
    throw InitException() << "RPVFFWVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " N: " << _theN << " U: " << _theU << " V: "
			  << _theV << Exception::abortnow;
  // SM interactions
  if(_interactions==0 || _interactions==1) {
    // particles for outgoing W-
    // quarks
    for(int ix=1;ix<6;ix+=2) {
      for(int iy=2;iy<7;iy+=2) {
	bool isOff = iy/2 != (ix+1)/2;
	if ( isOff && _diagonal )
	  continue;
	addToList(-ix, iy, -24);
      }
    }
    // leptons
    for(int ix=11;ix<17;ix+=2) {
      int inu = model->majoranaNeutrinos() ? (ix+23)/2 : ix+1;
      addToList(-ix, inu, -24);
    }
    // particles for outgoing W+
    // quarks
    for(int ix=2;ix<7;ix+=2) {
      for(int iy=1;iy<6;iy+=2) {
	bool isOff = ix/2 != (iy+1)/2;
	if ( isOff && _diagonal )
	  continue;
	addToList(-ix, iy, 24);
      }
    }
    // leptons
    for(int ix=11;ix<17;ix+=2) {
      int inu = model->majoranaNeutrinos() ? (ix+23)/2 : -ix-1;
      addToList(inu, ix, 24);
    }
  }
  // neutralino and chargino
  if(_interactions==0 || _interactions==2) {
    vector<long> neu(4);
    neu[0] = 1000022; neu[1] = 1000023;
    neu[2] = 1000025; neu[3] = 1000035;
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
    vector<long> cha(2);
    cha[0] = 1000024; cha[1] = 1000037;
    if(_theV->size().first==5) {
      cha.push_back(-11);
      cha.push_back(-13);
      cha.push_back(-15);
    }
    // sign == -1 outgoing W-, sign == +1 outgoing W+
    for(int sign = -1; sign < 2; sign += 2) {
      for(unsigned int ine = 0; ine < neu.size(); ++ine) {
	for(unsigned int ic = 0; ic < cha.size(); ++ic ) {
	  if(ic>1&&ine>3&&ic==ine-2) continue;
	  addToList(-sign*cha[ic], neu[ine], sign*24);
	}
      }
    }
  }
  Helicity::FFVVertex::doinit();
  // CKM matric
  if ( !_diagonal ) {
    Ptr<CKMBase>::transient_pointer CKM = model->CKM();
    // cast the CKM object to the HERWIG one
    ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
      hwCKM = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
					transient_const_pointer>(CKM);
    if(hwCKM) {
      vector< vector<Complex > > CKM;
      CKM = hwCKM->getUnsquaredMatrix(generator()->standardModel()->families());
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  _ckm[ix][iy]=CKM[ix][iy];
	}
      }
    }
    else {
      throw Exception() << "Must have access to the Herwig::StandardCKM object"
			<< "for the CKM matrix in RPVFFWVertex::doinit()"
			<< Exception::runerror;
    }
  }
  _sw = sqrt(sin2ThetaW());
}

void RPVFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _theN << _theU << _theV  << _diagonal << _ckm << _interactions;
}

void RPVFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _theN >> _theU >> _theV  >> _diagonal >> _ckm >> _interactions;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<RPVFFWVertex,Helicity::FFVVertex>
  describeHerwigRPVFFWVertex("Herwig::RPVFFWVertex", "HwSusy.so HwRPV.so");

void RPVFFWVertex::Init() {

  static ClassDocumentation<RPVFFWVertex> documentation
    ("The couplings of the fermions to the W boson in the RPV model"
     " with bilinear R-parity violation");

  static Switch<RPVFFWVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVFFWVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include all the interactions",
     0);
  static SwitchOption interfaceInteractionsSM
    (interfaceInteractions,
     "SM",
     "Only include the MS terms",
     1);
  static SwitchOption interfaceInteractionsSUSY
    (interfaceInteractions,
     "SUSY",
     "Include the neutralino/chargino terms",
     2);

  static Switch<RPVFFWVertex,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &RPVFFWVertex::_diagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);

}

void RPVFFWVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			       tcPDPtr part2,
#ifndef NDEBUG
  tcPDPtr part3) {
#else
  tcPDPtr) {
#endif
  assert(abs(part3->id()) == ParticleID::Wplus);
  // normalization
  // first the overall normalisation
  if(q2 != _q2last||_couplast==0.) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // left and right couplings for quarks
  if(abs(part1->id()) <= 6) {
    int iferm=abs(part1->id());
    int ianti=abs(part2->id());
    if(iferm%2!=0) swap(iferm,ianti);
    iferm = iferm/2;
    ianti = (ianti+1)/2;
    assert( iferm>=1 && iferm<=3 && ianti>=1 && ianti<=3);
    left(-sqrt(0.5)*_ckm[iferm-1][ianti-1]);
    right(0.);
  }
  else {
    long neu, cha;
    if(part1->charged()) {
      cha = part1->id();
      neu = part2->id();
    }
    else {
      cha = part2->id();
      neu = part1->id();
    }
    if(_theV->size().first==2&&abs(neu)<=16) {
      left(-sqrt(0.5));
      right(0.);
    }
    else {
      if(cha != _id1last || neu != _id2last) {
        _id1last = cha;
        _id2last = neu;
	unsigned int eigc = RPV_helper::charginoIndex(cha);
	unsigned int eign = RPV_helper::neutralinoIndex(neu);
	_leftlast = (*_theN)(eign, 1)*conj((*_theV)(eigc, 0)) - 
	  ( (*_theN)(eign, 3)*conj((*_theV)(eigc, 1))/sqrt(2));
	_rightlast = conj((*_theN)(eign, 1))*(*_theU)(eigc, 0) +
	  ( conj((*_theN)(eign, 2))*(*_theU)(eigc, 1)/sqrt(2));
	if(_theV->size().first==5) {
	  for(unsigned int k=0;k<3;++k)
	    _rightlast += ( conj((*_theN)(eign, 4+k))*(*_theU)(eigc, 2+k)/sqrt(2));
	}
      }
      Complex ltemp = _leftlast;
      Complex rtemp = _rightlast;
      bool chapart = abs(cha)>1000000 ? cha>0 : cha<0;
      // conjugate if +ve chargino
      if(chapart) {
	ltemp = conj(ltemp);
	rtemp = conj(rtemp);
      }
      if((part1->id()==cha&&chapart)||(part2->id()==cha&&!chapart)) {
	Complex temp = ltemp;
	ltemp  = -rtemp;
	rtemp = -temp;
      }
      left (ltemp);
      right(rtemp);
    }
  }
}
