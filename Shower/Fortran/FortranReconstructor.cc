// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranReconstructor class.
//

#include "FortranReconstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"

using namespace Herwig;

NoPIOClassDescription<FortranReconstructor> FortranReconstructor::initFortranReconstructor;
// Definition of the static class description member.

void FortranReconstructor::Init() {

  static ClassDocumentation<FortranReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering,",
      " including the kinematics reshuffling necessary to compensate for the recoil"
      "of the emissions." );

}

bool FortranReconstructor::reconstructHardJets(ShowerTreePtr,
					       map<tShowerProgenitorPtr,
					       pair<Energy,double> > pt ) const {
  return true;
  //throw Exception() << "FortranReconstructor::reconstructHardJets()"
  //		    << " not implemented yet" << Exception::runerror;
}

bool FortranReconstructor::reconstructDecayJets(ShowerTreePtr decay) const {
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator mit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mjt;
  vector<ShowerProgenitorPtr> ShowerHardJets;
  ShowerProgenitorPtr incoming;
  // extract the decaying particle
  for(mit=decay->incomingLines().begin();mit!=decay->incomingLines().end();++mit)
    incoming=(*mit).first;
  // extract the decay products
  for(mjt=decay->outgoingLines().begin();mjt!=decay->outgoingLines().end();++mjt)
    ShowerHardJets.push_back((*mjt).first);
  // loop over them
  for(unsigned int ix=0;ix<ShowerHardJets.size();++ix) {
    // did it radiate
    if(!ShowerHardJets[ix]->progenitor()->children().empty()) {
      // reconstruct the mass of the jet
      reconstructTimeLikeMass(ShowerHardJets[ix]->progenitor());
    }
  }

  throw Exception() << "FortranReconstructor::reconstructDecayJets()"
		    << " not implemented yet" << Exception::runerror;
}

void FortranReconstructor::reconstructTimeLikeMass(tShowerParticlePtr parent) const {
  // loop over the children and calculate the masses
  ParticleVector::const_iterator cit;
  for(cit=parent->children().begin();cit!=parent->children().end();++cit) {
    tShowerParticlePtr theLast=dynamic_ptr_cast<tShowerParticlePtr>(*cit);
    if((**cit).children().empty()) {
      parent->showerKinematics()->reconstructLast(theLast,0);
    }
    else {
      reconstructTimeLikeMass(theLast);
    }
  }
  // compute the mass of the particle
  parent->showerKinematics()->reconstructParent(parent,parent->children());
// C---SPECIAL FOR G-->QQBAR: READJUST ANGULAR DISTRIBUTION
//             IF (IDPAR(IPAR).EQ.13 .AND. IDPAR(JPAR).LT.13) THEN
//               Z=PPAR(4,JPAR)/PPAR(4,IPAR)
//               ZMIN=HWBVMC(IDPAR(JPAR))/PPAR(1,JPAR)*Z
//               RHO=(Z*(3-Z*(3-2*Z))-ZMIN*(3-ZMIN*(3-2*ZMIN)))
//      $             /(2*(1-2*ZMIN)*(1-ZMIN*(1-ZMIN)))
//               NQ=PPAR(3,IPAR)*(PPAR(3,IPAR)+PPAR(4,IPAR))
//               EMI=PPAR(5,IPAR)
//               EMJ=PPAR(5,JPAR)
//               EMK=PPAR(5,KPAR)
//               ZMIN=MAX((EMI+EMJ-EMK)/(2*(EMI+NQ)),
//      $      (EMI+EMJ-EMK-SQRT(ABS((EMI-EMJ-EMK)**2-4*EMJ*EMK)))/(2*EMI))
//               ZMAX=1-MAX((EMI-EMJ+EMK)/(2*(EMI+NQ)),
//      $      (EMI-EMJ+EMK-SQRT(ABS((EMI-EMJ-EMK)**2-4*EMJ*EMK)))/(2*EMI))
//               C=2*RMASS(IDPAR(JPAR))**2/EMI
//               Z=(4*ZMIN*(1.5*(1+C-ZMIN)+ZMIN**2)*(1-RHO)
//      $          +4*ZMAX*(1.5*(1+C-ZMAX)+ZMAX**2)*RHO-2-3*C)/(1+2*C)**1.5
//               Z=SQRT(1+2*C)*SINH(LOG(Z+SQRT(Z**2+1))/3)+0.5
//               Z=(Z*NQ+(EMI+EMJ-EMK)/2)/(NQ+EMI)
//               PPAR(4,JPAR)=Z*PPAR(4,IPAR)
//               PPAR(4,KPAR)=PPAR(4,IPAR)-PPAR(4,JPAR)
//               PPAR(3,JPAR)=HWUSQR(PPAR(4,JPAR)**2-EMJ)
//               PPAR(3,KPAR)=HWUSQR(PPAR(4,KPAR)**2-EMK)
//               PPAR(2,JPAR)=EXI/(PPAR(4,JPAR)*PPAR(4,KPAR))
//               IF(JDAPAR(2,JPAR).NE.0)PPAR(2,JDAPAR(2,JPAR))=PPAR(2,JPAR)
//               IF(JDAPAR(2,KPAR).NE.0)PPAR(2,JDAPAR(2,KPAR))=PPAR(2,JPAR)
// C---FIND DESCENDENTS OF THIS SPLITTING AND READJUST THEIR MOMENTA TOO
//               DO 20 J=JPAR+2,NPAR-1,2
//                 I=J
//  10             I=JMOPAR(1,I)
//                 IF (I.GT.IPAR) GOTO 10
//                 IF (I.EQ.IPAR) THEN
//                   I=JMOPAR(1,J)
//                   K=J+1
//                   POLD=PPAR(3,J)+PPAR(3,K)
//                   EOLD=PPAR(4,J)+PPAR(4,K)
//                   PNEW=HWUSQR(PPAR(4,I)**2-PPAR(5,I))
//                   ENEW=PPAR(4,I)
//                   A=(ENEW*EOLD-PNEW*POLD)/PPAR(5,I)
//                   B=(PNEW*EOLD-ENEW*POLD)/PPAR(5,I)
//                   PPAR(3,J)=A*PPAR(3,J)+B*PPAR(4,J)
//                   PPAR(4,J)=(PPAR(4,J)+B*PPAR(3,J))/A
//                   PPAR(3,K)=PNEW-PPAR(3,J)
//                   PPAR(4,K)=ENEW-PPAR(4,J)
//                   PPAR(2,J)=1-(PPAR(3,J)*PPAR(3,K)+PPAR(1,J)*PPAR(1,K))
//      $                 /(PPAR(4,J)*PPAR(4,K))
//                   IF (JDAPAR(2,J).NE.0) PPAR(2,JDAPAR(2,J))=PPAR(2,J)
//                   IF (JDAPAR(2,K).NE.0) PPAR(2,JDAPAR(2,K))=PPAR(2,J)
//                 ENDIF
//  20           CONTINUE
//             ENDIF


  
}

// CDECK  ID>, HWBMAS.
// *CMZ :-        -26/04/91  11.11.54  by  Bryan Webber
// *-- Author :    Bryan Webber
// C-----------------------------------------------------------------------
//       SUBROUTINE HWBMAS
// C-----------------------------------------------------------------------
// C     Passes  backwards through a  jet cascade  calculating the masses
// C     and magnitudes of the longitudinal and transverse three momenta.
// C     Components given relative to direction of parent for a time-like
// C     vertex and with respect to z-axis for space-like vertices.
// C
// C     On input PPAR(1-5,*) contains:
// C     (E*sqrt(Xi),Xi,3-mom (if external),E,M-sq (if external))
// C
// C     On output PPAR(1-5,*) (if TMPAR(*)), containts:
// C     (P-trans,Xi or Xilast,P-long,E,M)
// C-----------------------------------------------------------------------
//       INCLUDE 'HERWIG65.INC'
//       DOUBLE PRECISION HWUSQR,EXI,PISQ,PJPK,EJEK,PTSQ,Z,ZMIN,ZMAX,
//      $     EMI,EMJ,EMK,C,NQ,HWBVMC,RHO,POLD,PNEW,EOLD,ENEW,A,B
//       INTEGER IPAR,JPAR,KPAR,MPAR,I,J,K
//       EXTERNAL HWUSQR
//       IF (IERROR.NE.0) RETURN
//       IF (NPAR.GT.2) THEN
//         DO 30 MPAR=NPAR-1,3,-2
//          JPAR=MPAR
// C Find parent and partner of this branch
//          IPAR=JMOPAR(1,JPAR)
//          KPAR=JPAR+1
// C Determine type of branching
//          IF (TMPAR(IPAR)) THEN
//          ELSE
// C Space-like branching
// C           Re-arrange such that JPAR is time-like
//             IF (TMPAR(KPAR)) THEN
//                KPAR=JPAR
//                JPAR=JPAR+1
//             ENDIF
// C           Compute time-like branch
//             PTSQ=(2.-PPAR(2,JPAR))*PPAR(1,JPAR)*PPAR(1,JPAR)
//      &          -PPAR(5,JPAR)
//             PPAR(1,JPAR)=HWUSQR(PTSQ)
//             PPAR(3,JPAR)=(1.-PPAR(2,JPAR))*PPAR(4,JPAR)
//             PPAR(3,IPAR)=PPAR(3,KPAR)-PPAR(3,JPAR)
//             PPAR(5,IPAR)=0.
//             PPAR(1,KPAR)=0.
//          ENDIF
// C Reset Xi to Xilast
//          PPAR(2,KPAR)=PPAR(2,IPAR)
//  30    CONTINUE
//       ENDIF
//       DO 40 IPAR=2,NPAR
//  40   PPAR(5,IPAR)=HWUSQR(PPAR(5,IPAR))
//       PPAR(1,2)=0.
//       PPAR(2,2)=0.
//       END

bool FortranReconstructor::deconstructDecayJets(HardTreePtr,
						  EvolverPtr) const {
  throw Exception() << "FortranReconstructor::deconstructDecayJets() not "
		    << "implemented " << Exception::runerror;
}

bool FortranReconstructor::deconstructHardJets(HardTreePtr, EvolverPtr) const {
  throw Exception() << "FortranReconstructor::deconstructHardJets() "
		    << "not yet implemented"
		    << Exception::runerror;
}
