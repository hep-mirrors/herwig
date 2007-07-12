// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the UEDW0W1W1Vertex class.
//

#include "UEDW0W1W1Vertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

UEDW0W1W1Vertex::UEDW0W1W1Vertex() : theSinW(0.), theCosW(0.),
				     theSinThetaOne(0.), theCosThetaOne(0.),
				     theq2last(), theElast(0.), theCouplast(0.),
				     theSMlast(0), theKKlast(0) {
  vector<int> first(6), second(6), third(6);
  first[0] = 22;
  second[0] = -5100024;
  third[0] = 5100024;

  first[1] = 24;
  second[1] = -5100024;
  third[1] = 5100022;

  first[2] = -24;
  second[2] = 5100024;
  third[2] = 5100022;

  first[3] = 24;
  second[3] = -5100024;
  third[3] = 5100023;

  first[4] = -24;
  second[4] = 5100024;
  third[4] = 5100023;

  first[5] = 23;
  second[5] = -5100024;
  third[5] = 5100024;

  setList(first, second, third);
}

void UEDW0W1W1Vertex::doinit() throw(InitException) {
  VVVVertex::doinit();
   theUEDBase = dynamic_ptr_cast<tUEDBasePtr>(generator()->standardModel());
  if(!theUEDBase)
    throw InitException() << "UEDW0W1W1Vertex::doinit() - The pointer to "
			  << "the UEDBase object is null!"
			  << Exception::runerror;
  theSinW = sqrt(theUEDBase->sin2ThetaW());
  theCosW = sqrt( 1. - sqr(theSinW) );
  theSinThetaOne = theUEDBase->sinThetaOne();
  theCosThetaOne = sqrt( 1. - sqr(theSinThetaOne));
  orderInGs(0);
  orderInGem(1);
}

void UEDW0W1W1Vertex::persistentOutput(PersistentOStream & os) const {
  os << theUEDBase << theSinW << theCosW << theSinThetaOne 
     << theCosThetaOne;
}

void UEDW0W1W1Vertex::persistentInput(PersistentIStream & is, int) {
  is >> theUEDBase >> theSinW >> theCosW >> theSinThetaOne 
     >> theCosThetaOne;
}

ClassDescription<UEDW0W1W1Vertex> UEDW0W1W1Vertex::initUEDW0W1W1Vertex;
// Definition of the static class description member.

void UEDW0W1W1Vertex::Init() {

  static ClassDocumentation<UEDW0W1W1Vertex> documentation
    ("The coupling of an SM W boson to a level 1 KK W and KK Z and KK photon");

}

void UEDW0W1W1Vertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				  tcPDPtr part3) {
  long id1(abs(part1->id())), id2(abs(part2->id())), id3(abs(part3->id())), 
    smID(0), kkparticle(0);
  double perm(-1.);
  if( id1 == 22 || id1 == 23 || id1 == 24 ) {
    if( part1->id() < 0 ) perm = 1.;
    smID = id1;
    kkparticle = (id2 == 5100024) ? id3 : id2;
  }
  else if( id2 == 22 || id2 == 23 || id2 == 24 ) {
    if( part2->id() < 0 ) perm = 1.;
    smID = id2;
    kkparticle = (id1 == 5100024) ? id3 : id1;
  }
  else if( id3 == 22 || id3 == 23 || id3 == 24 ) {
    if( part3->id() < 0 ) perm = 1.;
    smID = id3;
    kkparticle = (id1 == 5100024) ? id2 : id1;
  }
  else {
    throw HelicityLogicalError()
      << "UEDW0W1W1Vertex::setCoupling() - There is no SM gauge boson in "
      << "this vertex. " << id1 << " " << id2 << " " << id3 
      << Exception::warning; 
    setNorm(0.);
    return;
  }
  if( q2 != theq2last ) {
    theq2last = q2;
    theElast = sqrt(4.*Constants::pi*theUEDBase->alphaEM(q2));
  }
  
  if( smID != theSMlast || kkparticle != theKKlast ) { 
    theSMlast = smID;
    theKKlast = kkparticle;
    if( smID == 22 )
      theCouplast = 1.;
    else if(smID == 23) 
      theCouplast = theCosW/theSinW;
    else {
      if( kkparticle == 5100023 )
	theCouplast = theCosThetaOne/theSinW;
      else
	theCouplast = theSinThetaOne/theSinW;
    }
  }
  setNorm(perm*theElast*theCouplast);
}
