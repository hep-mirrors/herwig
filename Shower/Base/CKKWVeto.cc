// -*- C++ -*-
//
// CKKWVeto.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CKKWVeto class.
//

#include "CKKWVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ShowerProgenitor.h"
#include "Herwig++/Shower/ShowerHandler.h"

using namespace Herwig;

void CKKWVeto::persistentOutput(PersistentOStream & os) const {
  os <<  _vetoTimelike << _vetoSpacelike << _ptVetoDefinition << ounit(_Ptveto,GeV) 
     << _reversePtVeto;
}

void CKKWVeto::persistentInput(PersistentIStream & is, int) {
  is  >> _vetoTimelike >> _vetoSpacelike >> _ptVetoDefinition >> iunit(_Ptveto,GeV) 
      >> _reversePtVeto;
}

ClassDescription<CKKWVeto> CKKWVeto::initCKKWVeto;
// Definition of the static class description member.

void CKKWVeto::Init() {

  static ClassDocumentation<CKKWVeto> documentation
    ("CKKW vetoing in shower");
   
  static Switch<CKKWVeto, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &CKKWVeto::_ptVetoDefinition, 1, false, false);
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet measure", 0);
  static SwitchOption Shower
    (ifaceJetMeasureMode,"Shower","Shower pt", 1);
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet measure", 2);
  static SwitchOption Hadron
    (ifaceJetMeasureMode,"Hadron","Hadron jet measure", 3);
 

  static Parameter< CKKWVeto, Energy > interfacePtCut
    ("JetCut",
     "The jet cut (in specified jet measure) for shower vetoes",
     &CKKWVeto::_Ptveto, GeV, ZERO, ZERO, 100000.0 * GeV,
     false, false, Interface::limited );

  static Switch<CKKWVeto, bool> ifaceReverseVeto
    ("ReversePtVeto",
     "Reverse pt veto to veto emissions below cut",
     &CKKWVeto::_reversePtVeto, false, false, false);
  static SwitchOption RevVetoFalse
    (ifaceReverseVeto,"No","Veto emissions above cut", false);
  static SwitchOption RevVetoTrue
    (ifaceReverseVeto,"Yes","Veto emissions below cut", true);

  static Switch<CKKWVeto,bool> interfacevetoTimelike
    ("VetoTimelike",
     "Veto timelike showering",
     &CKKWVeto::_vetoTimelike, true, false, false);
  static SwitchOption interfacevetoTimelikevetoTimelikeOn
    (interfacevetoTimelike,
     "Yes",
     "Veto timelike showering",
     true);
  static SwitchOption interfacevetoTimelikevetoTimelikeOff
    (interfacevetoTimelike,
     "No",
     "Do not veto timelike showering",
     false);

  static Switch<CKKWVeto,bool> interfacevetoSpacelike
    ("VetoSpacelike",
     "Veto spacelike showering",
     &CKKWVeto::_vetoSpacelike, true, false, false);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOn
    (interfacevetoSpacelike,
     "Yes",
     "Veto spacelike showering",
     true);
  static SwitchOption interfacevetoSpacelikevetoSpacelikeOff
    (interfacevetoSpacelike,
     "No",
     "Do not veto spacelike showering",
     false);
}

bool CKKWVeto::vetoTimeLike (tcShowerProgenitorPtr progenitor, tcShowerParticlePtr,
			     const Branching & fb ){
  //find pt at which we are vetoing
  Energy2 ptVeto;
  if( _Ptveto > ZERO && !_highestMult ) ptVeto = sqr( _Ptveto );
  else ptVeto = sqr( progenitor->maximumpT() );
  //find corresponding pt measure based on the emission variables
  Energy2 kt_measure;
  Energy2 s = ShowerHandler::currentHandler()->lastXCombPtr()->lastS();
  // Energy2 s = sqr( 91.2*GeV);
  Energy pt = fb.kinematics->pT();
  double z = fb.kinematics->z();
  //exact single durham/luclus cuts
  if( fb.kinematics && ( _ptVetoDefinition == 0 || _ptVetoDefinition == 2 ) 
      && ! _highestMult ){
    Energy2 m0 = sqr(getParticleData( fb.ids[0] )->constituentMass());
    Energy2 m1 = sqr(getParticleData( fb.ids[1] )->constituentMass());
    Energy2 m2 = sqr(getParticleData( fb.ids[2] )->constituentMass());
    
    double lambda = sqrt( 1. - 4.*m0/s );
    double beta1 = 2.*( m1 - sqr(z)*m0 + sqr(pt) )
      / z / lambda / ( lambda + 1. ) / s;
    double beta2 = 2.*( m2 - sqr( 1. - z )*m0 + sqr(pt) )
      / ( 1. - z ) / lambda / ( lambda + 1. ) / s;
    
    Energy E1 = sqrt(s)/2.*( z + lambda*beta1 );
    Energy E2 = sqrt(s)/2.*( (1.-z) + lambda*beta2 );
    Energy Z1 = sqrt(s)/2.*lambda*( z - beta1 );
    Energy Z2 = sqrt(s)/2.*lambda*( (1.-z) - beta2 );
    
    double costheta = ( Z1*Z2 - sqr(pt) )
      / sqrt( sqr(Z1)+sqr(pt) ) / sqrt( sqr(Z2)+sqr(pt) );
    
    if( _ptVetoDefinition == 0 )
      kt_measure = 2.*min( sqr(E1), sqr(E2) )*( 1. - costheta );
    else if( _ptVetoDefinition == 2 )
      kt_measure = 2.*sqr(E1)*sqr(E2)/sqr(E1+E2)*( 1. - costheta );
  }
  //hadron jet measure cuts
  else if( fb.kinematics && _ptVetoDefinition == 3 && !_highestMult ){
    Energy2 m1 = sqr(getParticleData( fb.ids[1] )->constituentMass());
    Energy2 m2 = sqr(getParticleData( fb.ids[2] )->constituentMass());
    
    double beta1 = 2.*( m1 + sqr(pt) ) / z  / s;
    double beta2 = 2.*( m2 + sqr(pt) ) / ( 1. - z ) / s;
    
    Energy E1 = sqrt(s)/2.*( z + beta1 );
    Energy E2 = sqrt(s)/2.*( (1.-z) + beta2 );
    Energy Z1 = sqrt(s)/2.*( z - beta1 );
    Energy Z2 = sqrt(s)/2.*( (1.-z) - beta2 );
      
    //delta phi is always pi for first emission (qt_i = +-pt)
    double deltaR = sqr(  log( z / beta1 ) - log( (1-z) / beta2 ) ) / 4. 
      + sqr( Constants::pi );
    kt_measure = sqr( pt )* deltaR;
  } 
  //normal shower pt veto - should always be called for highest mult
  else kt_measure = sqr( fb.kinematics->pT() );
  
  //veto emission or the full event 
  if( _dynamicSuds ){
    if( ! _reversePtVeto && kt_measure > ptVeto ) throw Veto();
    else if(  _reversePtVeto && kt_measure < ptVeto ) throw Veto();
  }
  else{
    if( ! _reversePtVeto && kt_measure > ptVeto ) return true;
    else if(  _reversePtVeto && kt_measure < ptVeto ) return true;
  }
  return false;
}

/**
 * Return true, if the selected emission off the given
 * particle and progenitor is vetoed.
 */
bool CKKWVeto::vetoSpaceLike (tcShowerProgenitorPtr progenitor, tcShowerParticlePtr,
			      const Branching & bb){
  Energy ptVeto;
  if( _Ptveto > ZERO && !_highestMult ) ptVeto = _Ptveto;
  else  ptVeto = progenitor->maximumpT();

  Energy kt_measure = bb.kinematics->pT();

  //veto emission or the full event 
  if( _dynamicSuds ){
    if( ! _reversePtVeto && kt_measure > ptVeto ) throw Veto();
    else if(  _reversePtVeto && kt_measure < ptVeto ) throw Veto();
  }
  else{
    if( ! _reversePtVeto && kt_measure > ptVeto ) return true;
    else if(  _reversePtVeto && kt_measure < ptVeto ) return true;
  }
  return false;
}
