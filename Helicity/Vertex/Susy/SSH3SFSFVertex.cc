// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSH3SFSFVertex class.
//

#include "SSH3SFSFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

SSH3SFSFVertex::SSH3SFSFVertex():_mw(0.),_sw(0.),_tb(0.),_mu(0.),
				 _trilinear(9,0.),_q2last(0.),
				 _glast(0.),_hfact(0.),_id1last(0),
				 _id2last(0) {
  vector<int> first,second,third;
  for(unsigned int ix=1000001;ix<1000007;++ix){
    first.push_back(36);
    second.push_back(ix);
    third.push_back(-(ix+1000000));
  }
  for(unsigned int ix=1000001;ix<1000007;++ix){
    first.push_back(36);
    second.push_back(-ix);
    third.push_back(ix+1000000);
  }
  for(unsigned int ix=1000011;ix<1000016;ix+=2){
    first.push_back(36);
    second.push_back(ix);
    third.push_back(-(ix+1000000));
  }
  for(unsigned int ix=1000011;ix<1000016;ix+=2){
    first.push_back(36);
    second.push_back(-ix);
    third.push_back(ix+1000000);
  }
  setList(first,second,third);
}

void SSH3SFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _mw << _sw << _tb << _mu << _trilinear;
}

void SSH3SFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _mw >> _sw >> _tb >> _mu >> _trilinear;
  _q2last=0.;
  _glast=0.;
  _hfact=0.;
  _id1last=0;
  _id2last=0;
}

void SSH3SFSFVertex::doinit() throw(InitException) {
  SSSVertex::doinit();
  _theSS = dynamic_ptr_cast<SusyBasePtr>(generator()->standardModel());
  if(!_theSS)
    throw InitException() << "SSGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _mw = getParticleData(24)->mass();
  _sw = _theSS->sin2ThetaW();
  _tb = _theSS->tanBeta();
  _mu = _theSS->muParameter();
  vector<Complex> au(_theSS->upTypeTrilinear());
  vector<Complex> ad(_theSS->downTypeTrilinear());
  vector<Complex> ae(_theSS->leptonTypeTrilinear());
  _trilinear[0] = ad[0].real();
  _trilinear[1] = au[0].real();
  _trilinear[2] = ad[1].real();
  _trilinear[3] = au[1].real();
  _trilinear[4] = ad[2].real();
  _trilinear[5] = au[2].real();
  _trilinear[6] = ae[0].real();
  _trilinear[7] = ae[1].real();
  _trilinear[8] = ae[2].real();
  orderInGem(1);
  orderInGs(0);
  }



ClassDescription<SSH3SFSFVertex> SSH3SFSFVertex::initSSH3SFSFVertex;
// Definition of the static class description member.

void SSH3SFSFVertex::Init() {

  static ClassDocumentation<SSH3SFSFVertex> documentation
    ("This class implements the coupling of the H_3 neutral higgs to "
     "sfermion pairs.");

}

 void SSH3SFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				  tcPDPtr part2,tcPDPtr part3) {
   long isf1(0),isf2(0);
   if(abs(part1->id()) == 36 ) {
     isf1 = abs(part2->id());
    isf2 = abs(part3->id());
   }
   else if(abs(part2->id()) == 36 ) {
     isf1 = abs(part1->id());
     isf2 = abs(part3->id());
   }
   else if(abs(part3->id()) == 36 ) {
     isf1 = abs(part1->id());
     isf2 = abs(part2->id());
   }
   else {
     throw HelicityConsistencyError() << "SSH2SSVertex::setCoupling "
				      << "No A0 particle found!\n"
				      << Exception::warning;
     setNorm(0.);
}
   
   if((isf1>=1000001 && isf1<=1000006)||(isf1>=2000001 && isf1<=2000006)||
      (isf1>=1000011 && isf1<=1000016)||(isf1>=2000011 && isf1<=2000016) &&
      isf1!=isf2) {
     if(q2 != _q2last) {
       double alpha = _theSS->alphaEM(q2);
       _glast = -Complex(1.,0.)*sqrt(4.*Constants::pi*alpha)/_sw;
       _q2last = q2;
     }
     if(isf1!=_id1last ||isf2!=_id2last) {
       _id1last = isf1;
       _id2last = isf2;
       unsigned int imass(0);
       if((isf1/1000000)==1) {imass=isf1-1000000;}
       else {imass=isf1-2000000;}
       Energy fmass = getParticleData(imass)->mass();
       _hfact = fmass/2./_mw;
       if(isf1%2==0) {
	 unsigned int tri;
	 if((isf1/1000000)==1) {tri=isf1-1000001;}
	 else {tri=isf1-2000001;}
	 _hfact *= (_mu + (_trilinear[tri]/_tb));
       }
       else {
	 unsigned int tri;
	 if(isf1>=1000011||isf1>=2000011) {
	   if((isf1/1000000==1)) {tri=(isf1-999999)/2;}
	   else{tri=(isf1-1999999)/2;}
	 }
	 else {
	   if((isf1/1000000)==1) {tri=isf1-1000001;}
	   else {tri=isf1-2000001;}
	 }
	 _hfact *= (_mu + _trilinear[tri]*_tb);
       }
     }
      setNorm(_glast*_hfact);
   }
   else {
     throw HelicityConsistencyError() << "Incorrect particles detected in "
				       << "SSH3SFSFVertex. "
				       << Exception::warning;
   }
 }
