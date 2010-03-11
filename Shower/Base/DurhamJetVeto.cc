// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DurhamJetVeto class.
//

#include "DurhamJetVeto.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr DurhamJetVeto::clone() const {
  return new_ptr(*this);
}

IBPtr DurhamJetVeto::fullclone() const {
  return new_ptr(*this);
}

void DurhamJetVeto::persistentOutput(PersistentOStream & os) const {
  os << _ptVetoDefinition << _reversePtVeto << _y_cut << _approxCuts;
}

void DurhamJetVeto::persistentInput(PersistentIStream & is, int) {
  is >> _ptVetoDefinition >> _reversePtVeto >> _y_cut >> _approxCuts;
}

ClassDescription<DurhamJetVeto> DurhamJetVeto::initDurhamJetVeto;
// Definition of the static class description member.

void DurhamJetVeto::Init() {

  static ClassDocumentation<DurhamJetVeto> documentation
    ("There is no documentation for the DurhamJetVeto class");
  
  static Switch<DurhamJetVeto, bool> ifaceApproxCuts
    ("ApproxCuts",
     "Use approximate rather than exact jet cuts",
     &DurhamJetVeto::_approxCuts, false, false, false);
  static SwitchOption ApproxCutsOff
    (ifaceApproxCuts,"No","Use exact cuts", false);
  static SwitchOption ApproxCutsOn
    (ifaceApproxCuts,"Yes","Use approx cuts", true);

  static Switch<DurhamJetVeto, bool> ifaceReverseVeto
    ("ReversePtVeto",
     "Reverse pt veto to veto emissions below cut",
     &DurhamJetVeto::_reversePtVeto, false, false, false);
  static SwitchOption RevVetoFalse
    (ifaceReverseVeto,"No","Veto emissions above cut", false);
  static SwitchOption RevVetoTrue
    (ifaceReverseVeto,"Yes","Veto emissions below cut", true);

  static Parameter<DurhamJetVeto,double> interfacePtCut
    ("JetCut",
     "The jet cut",
     &DurhamJetVeto::_y_cut, 1.1, 0., 1.1,
     false, false, Interface::limited );

  static Switch<DurhamJetVeto, unsigned int> ifaceJetMeasureMode
    ("JetMeasure",
     "Choice of the jet measure algorithm",
     &DurhamJetVeto::_ptVetoDefinition, 1, false, false);
  static SwitchOption Durham
    (ifaceJetMeasureMode,"Durham","Durham jet algorithm", 0);
  static SwitchOption Shower
    (ifaceJetMeasureMode,"Shower","Shower pt", 1);
  static SwitchOption LUCLUS
    (ifaceJetMeasureMode,"LUCLUS","LUCLUS jet algorithm", 2);

}

bool DurhamJetVeto::vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
		   const Branching & fb) {
  Energy2 s = ShowerHandler::currentHandler()->lastXCombPtr()->lastS();
  Energy ptVeto = sqrt( s * _y_cut );
  if( _ptVetoDefinition != 1 ) {
    Energy2 kt_measure(ZERO);
    double z = fb.kinematics->z();
    Energy pt = fb.kinematics->pT();
    if( !_approxCuts ) {
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
    else {
      
      if( _ptVetoDefinition == 0 )
	kt_measure = sqr( pt / max( z, 1. - z ) );
      else if( _ptVetoDefinition == 2 )
	kt_measure = sqr( pt );

    }
    if( ! _reversePtVeto ){ 
      if( kt_measure > sqr(ptVeto) ) return true;
    }
    else{
      if( kt_measure < sqr(ptVeto) ) return true;
    }
  }
  else {
    if(fb.kinematics->pT()> ptVeto) return true;
  }
  return false;
}

bool DurhamJetVeto::vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
		    const Branching& bb) {
  Energy ptVeto = sqrt( ShowerHandler::currentHandler()->lastXCombPtr()->lastS() 
			* _y_cut );
  if( bb.kinematics && _ptVetoDefinition == 0 ) 
    ptVeto *= max( bb.kinematics->z(), 1. - bb.kinematics->z() );
  if(bb.kinematics->pT()> ptVeto )
    return true;
  return false;
}
