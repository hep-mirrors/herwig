// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWSSVertex class.
//

#include "SSWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

inline SSWSSVertex::SSWSSVertex():_sw(0.),_q2last(),_couplast(0.) {
  vector<int> first,second,third;
  //W-
  //LL-squarks
  for(unsigned int ix=1000001;ix<1000006;ix+=2) {
    first.push_back(-24);
    second.push_back(ix+1);
    third.push_back(-ix);
  }
  //1-2 stop sbottom
  first.push_back(-24);
  second.push_back(1000006);
  third.push_back(-2000005);
  //2-1 stop sbottom
  first.push_back(-24);
  second.push_back(2000006);
  third.push_back(-1000005);
  //2-2 stop sbottom
  first.push_back(-24);
  second.push_back(2000006);
  third.push_back(-2000005);
 
  //LL-sleptons
  for(unsigned int ix=1000011;ix<1000016;ix+=2) {
    first.push_back(-24);
    second.push_back(-ix);
    third.push_back(ix+1);
  }
  //2-L stau
  first.push_back(-24);
  second.push_back(-2000015);
  third.push_back(1000016);
  //W+
  for(unsigned int ix=1000001;ix<1000006;ix+=2) {
    first.push_back(24);
    second.push_back(-(ix+1));
    third.push_back(ix);
  }

//1-2 stop sbottom
  first.push_back(24);
  second.push_back(-1000006);
  third.push_back(2000005);
  //2-1 stop sbottom
  first.push_back(24);
  second.push_back(-2000006);
  third.push_back(1000005);
  //2-2 stop sbottom
  first.push_back(24);
  second.push_back(-2000006);
  third.push_back(2000005);

  //LL-sleptons
  for(unsigned int ix=1000011;ix<1000016;ix+=2) {
    first.push_back(24);
    second.push_back(ix);
    third.push_back(-ix-1);
  }
  //2-L stau
  first.push_back(24);
  second.push_back(2000015);
  third.push_back(-1000016);
  
  setList(first,second,third);
}

void SSWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sw << _stau << _stop << _sbottom;
}

void SSWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _sw >> _stau >> _stop >> _sbottom;
  _q2last=0.*GeV2;
  _couplast=0.;
}

ClassDescription<SSWSSVertex> SSWSSVertex::initSSWSSVertex;
// Definition of the static class description member.

void SSWSSVertex::Init() {

  static ClassDocumentation<SSWSSVertex> documentation
    ("This is the implementation of the coupling of the W to two "
     "sfermions");
  
}

void SSWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3){
  long utype(0),dtype(0);//utype for leptons means sneutrinos
  if(abs(part1->id()) == 24) {
    if(part2->id() % 2 == 0) {
      utype = abs(part2->id());
      dtype = abs(part3->id());
    }
    else {
      utype = abs(part3->id());
      dtype = abs(part2->id());
    }
  }
  else if(abs(part2->id()) == 24) {
    if(part1->id() % 2 == 0) {
      utype = abs(part1->id());
      dtype = abs(part3->id());
    }
    else {
      utype = abs(part3->id());
      dtype = abs(part1->id());
    }
  }
  else {
    if(part1->id() % 2 == 0) {
      utype = abs(part1->id());
      dtype = abs(part2->id());
    }
    else {
      utype = abs(part2->id());
      dtype = abs(part1->id());
    }
  }
  unsigned int eig1(utype/1000000 -1),eig2(dtype/1000000 -1);
  if((eig1 == eig2) || utype==1000006 || utype==200006 ||
     dtype==1000015 || dtype==2000015) {
    if(q2 != _q2last) {
      _q2last = q2;
      double alpha = _theSS->alphaEM(q2); 
      _couplast = sqrt(2.*Constants::pi*alpha)/_sw;
      
      if(utype==1000006 ||utype==2000006){
	_couplast *= (*_stop)(0,eig1)*(*_sbottom)(0,eig2); 
      }
      if(dtype==1000015||dtype==2000015) {
	_couplast *= (*_stau)(0,eig2);
      }
    }
    
    setNorm(_couplast);
  }
  else {
    throw HelicityConsistencyError() << "SSWSSVertex::setCoupling "
				     << "Unknown particle "
				     << part1->PDGName() << " " 
				     << part2->PDGName() << " "
				     << part3->PDGName() 
				     << " in WSS vertex\n"
				     << Exception::warning;
    setNorm(0.);
  }
}

