// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSPSSVertex class.
//

#include "SSPSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSPSSVertex::SSPSSVertex():_couplast(0.),_q2last() {
  vector<int> first,second,third;
  //sleptons
  for(unsigned int ix=1000011;ix<1000016;ix+=2) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=2000011;ix<2000016;ix+=2) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //squarks
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=2000001;ix<2000007;++ix) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  setList(first,second,third);
}

void SSPSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS;
}

void SSPSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS;
  _couplast=0.;
  _q2last=0*GeV2;
}

ClassDescription<SSPSSVertex> SSPSSVertex::initSSPSSVertex;
// Definition of the static class description member.

void SSPSSVertex::Init() {

  static ClassDocumentation<SSPSSVertex> documentation
    ("The SSPSSVertex class implements the coupling of a photon "
     "to 2 sfermions.");
  
}

void SSPSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2, tcPDPtr) {
  long isf(0);
  if(part1->id()==22) {
    isf = abs(part2->id());
  }
  else if(part2->id()==22) {
    isf = abs(part1->id());
  }
  else {
    isf = abs(part1->id());
  }
  
  if((isf>=1000001 && isf<=1000006)||(isf>=2000001 && isf<=2000006)||
     isf==1000011||isf==1000013||isf==1000015||isf==2000011||
     isf==2000013||isf==2000015){
    if(q2 != _q2last) {
      double alpha = _theSS->alphaEM(q2);
      _couplast = sqrt(4.*Constants::pi*alpha);
      if(isf>=1000011){
	_couplast *= -1.;
      }
      else if((isf<=1000006 && (isf%2)==0)||(isf<=2000006 && (isf%2)==0)){
	_couplast *= _theSS->eu();
      }
      else {
	_couplast *= _theSS->ed();
      }
      _q2last = q2;
    }
    setNorm(_couplast);

  }
  else {
    throw  HelicityConsistencyError() << "SSPSSVertex::setCoupling() - "
				      << "Incorrect particle(s) in " 
				      << "SSPSSVertex."
				      << Exception::warning;
  }

}
 


