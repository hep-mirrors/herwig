// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFWVertex class.
//

#include "SMFFWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
SMFFWVertex::SMFFWVertex() : _ckm(3,vector<Complex>(3,0.0)),_couplast(0.),
			     _q2last(0.*sqr(MeV))
{
  // particles for the vertex
  vector<int> first,second,third;
  // particles for outgoing W-
  // quarks
  for(unsigned int ix=1;ix<6;ix+=2) {
    for(unsigned int iy=2;iy<7;iy+=2) {
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(-24);
    }
  }
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix);
    second.push_back(ix+1);
    third.push_back(-24);
  }
  // particles for outgoing W+
  // quarks
  for(unsigned int ix=2;ix<7;ix+=2) {
    for(unsigned int iy=1;iy<6;iy+=2) {
      first.push_back(-ix);
      second.push_back(iy);
      third.push_back(24);
    }
  }
  // leptons
  for(unsigned int ix=11;ix<17;ix+=2) {
    first.push_back(-ix-1);
    second.push_back(ix);
    third.push_back(24);
  }
  setList(first,second,third);
}

void SMFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _theCKM << _theSM << _ckm;
}
  
void SMFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theCKM >> _theSM >> _ckm;
}
  
void SMFFWVertex::doinit() throw(InitException) {
  Herwig::Helicity::FFVVertex::doinit();
  _theSM  = generator()->standardModel();
  _theCKM = generator()->standardModel()->CKM();
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer>(_theCKM);
  if(hwCKM) {
    vector< vector<Complex > > CKM;
    CKM = hwCKM->getUnsquaredMatrix(_theSM->families());
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int iy=0;iy<3;++iy) {
	_ckm[ix][iy]=CKM[ix][iy];
      }
    }
  }
  else {
    throw InitException() << "Must have access to the Herwig++::StandardCKM object"
			  << "for the CKM matrix in SMFFWVertex::doinit()"
			  << Exception::runerror;
  }
  orderInGem(1);
  orderInGs(0);
}

ClassDescription<SMFFWVertex>SMFFWVertex::initSMFFWVertex;
  
// Definition of the static class description member.
  
void SMFFWVertex::Init() {
  static ClassDocumentation<SMFFWVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the W boson to the Standard Model fermions");
  
}
  
// coupling for FFW vertex
void SMFFWVertex::setCoupling(Energy2 q2, tcPDPtr a, tcPDPtr b, tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    double sw    = sqrt(2.*(_theSM->sin2ThetaW()));
    _couplast = -sqrt(4.0*3.14159265*alpha)/sw;
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  // quarks
  if(iferm>=1 && iferm <=6) {
    int iu,id;
    // up type first
    if(iferm%2==0) {
      iu = iferm/2;
      id = (ianti+1)/2;
    }
    // down type first
    else {
      iu = ianti/2;
      id = (iferm+1)/2;
    }
    if( iu<1 || iu>3 || id<1 || id>3) {
      throw HelicityConsistencyError() << "SMFFWVertex::setCoupling "
				       << "Unknown particle in W vertex" 
				       << Exception::warning;
      setLeft(0.);
      setRight(0.);
      return;
    }
    setLeft(_ckm[iu-1][id-1]);
    setRight(0.);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    setLeft(1.);
    setRight(0.);
  }
  else {
    throw HelicityConsistencyError() << "SMFFWVertex::setCoupling "
				     << "Unknown particle in W vertex" 
				     << Exception::warning;
    setLeft(0.);setRight(0.);
  }
}
 
}
}





