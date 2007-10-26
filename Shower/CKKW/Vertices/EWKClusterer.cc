// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EWKClusterer class.
//

#include "EWKClusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EWKClusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

EWKClusterer::~EWKClusterer() {}

void EWKClusterer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _useColour << _useHadronic << _EWKOrder;
}

void EWKClusterer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _useColour >> _useHadronic >> _EWKOrder;  
}

ClassDescription<EWKClusterer> EWKClusterer::initEWKClusterer;
// Definition of the static class description member.

void EWKClusterer::Init() {

  static ClassDocumentation<EWKClusterer> documentation
    ("Electroweak vertices in ME/PS merging");


  static Switch<EWKClusterer,bool> interfaceUseColour
    ("UseColour",
     "Wether or not to use colour information on clusterings",
     &EWKClusterer::_useColour, true, false, false);
  static SwitchOption interfaceUseColourUseColourOn
    (interfaceUseColour,
     "Yes",
     "Use colour on clustering",
     true);
  static SwitchOption interfaceUseColourUseColourOff
    (interfaceUseColour,
     "No",
     "Do not use colour on clustering",
     false);

  static Switch<EWKClusterer,bool> interfaceUseHadronic
    ("UseHadronic",
     "Wether or not to use hadronic clusterings",
     &EWKClusterer::_useHadronic, true, false, false);
  static SwitchOption interfaceUseHadronicUseHadronicOn
    (interfaceUseHadronic,
     "Yes",
     "Use hadronic clusterings",
     true);
  static SwitchOption interfaceUseHadronicUseHadronicOff
    (interfaceUseHadronic,
     "No",
     "Do not use hadronic clusterings",
     false);



  static Parameter<EWKClusterer,unsigned int> interfaceEWKOrder
    ("EWKOrder",
     "The electroweak order to consider",
     &EWKClusterer::_EWKOrder, 2, 0, 0,
     false, false, Interface::lowerlim);


}

ClusteringParticleData EWKClusterer::doEmergingLine
(ClusteringParticleData p1, ClusteringParticleData p2, bool& isEWK) const {

  if (!_useColour) throw Exception () << "CKKW : EWKClusterer::doEmergingLine : have to use colour information"
				      << Exception::runerror;

  ClusteringParticleData emerging;

  // only consider final state particles here

  if (p1.partonId.state == ClusteringParticleState::initial ||
      p2.partonId.state == ClusteringParticleState::initial) {
    isEWK = false;
    return emerging;
  }

  emerging.partonId.state = ClusteringParticleState::final;

  // quarks

  if (_useHadronic) {

    // cluster quark and Z

    if (abs(p1.partonId.PDGId) < 7 && p2.partonId.PDGId == 23) {
      emerging.partonId.PDGId = p1.partonId.PDGId;
      emerging.colour = p1.colour;
      emerging.antiColour = p1.antiColour;
      isEWK = true;
      return emerging;
    }
    
    // cluster quark and W+
    
    if (abs(p1.partonId.PDGId) < 7 && p2.partonId.PDGId == 24 &&
      ((p1.partonId.PDGId > 0 && p1.partonId.PDGId % 2 == 1) || // d, s, b
       (p1.partonId.PDGId < 0 && p1.partonId.PDGId % 2 == 0))) { // ubar, cbar, tbar
      emerging.partonId.PDGId = p1.partonId.PDGId + 1;
      emerging.colour = p1.colour;
      emerging.antiColour = p1.antiColour;
      isEWK = true;
      return emerging;
    }
    
    // cluster quark and W-
    
    if (abs(p1.partonId.PDGId) < 7 && p2.partonId.PDGId == -24 &&
	((p1.partonId.PDGId < 0 && p1.partonId.PDGId % 2 == 1) || // dbar, sbar, bbar
	 (p1.partonId.PDGId > 0 && p1.partonId.PDGId % 2 == 0))) { // u, c, t
      emerging.partonId.PDGId = p1.partonId.PDGId - 1;
      emerging.colour = p1.colour;
      emerging.antiColour = p1.antiColour;
      isEWK = true;
      return emerging;
    }

    // cluster quarks to Z
    
    if (abs(p1.partonId.PDGId) < 7 && p1.partonId.PDGId+p2.partonId.PDGId == 0 &&
      ((p1.colour == p2.antiColour) || (p1.antiColour == p2.colour))) {
      emerging.partonId.PDGId = 23;
      emerging.colour = 0;
      emerging.antiColour = 0;
      isEWK = true;
      return emerging;
    }

    // cluster quarks to W+

    if (abs(p1.partonId.PDGId)<7 && 
	p1.partonId.PDGId < 0 && p2.partonId.PDGId > 0 &&
	p1.partonId.PDGId + p2.partonId.PDGId == 1) {
      emerging.partonId.PDGId = 24;
      emerging.colour = 0;
      emerging.antiColour = 0;
      isEWK = true;
      return emerging;
    }

    // cluster quarks to W-
    
    if (abs(p1.partonId.PDGId)<7 && 
	p1.partonId.PDGId < 0 && p2.partonId.PDGId > 0 &&
	p1.partonId.PDGId + p2.partonId.PDGId == -1) {
      emerging.partonId.PDGId = -24;
      emerging.colour = 0;
      emerging.antiColour = 0;
      isEWK = true;
      return emerging;
    }

  } // if (_useHadronic)

  // only consider decay to leptons

  // l lbar to z

  if (abs(p1.partonId.PDGId)>10 && abs(p1.partonId.PDGId)<12 &&
      p1.partonId.PDGId+p2.partonId.PDGId == 0) {
    emerging.partonId.PDGId = 23;
    emerging.colour = 0;
    emerging.antiColour = 0;
    isEWK = true;
    return emerging;
  }

  // lbar nu to W+

  if (abs(p1.partonId.PDGId)>10 && abs(p1.partonId.PDGId)<12 &&
      p1.partonId.PDGId + p2.partonId.PDGId == 1) {
    emerging.partonId.PDGId = 24;
    emerging.colour = 0;
    emerging.antiColour = 0;
    isEWK = true;
    return emerging;
  }

  // l nubar to W-

  if (abs(p1.partonId.PDGId)>10 && abs(p1.partonId.PDGId)<12 &&
      p1.partonId.PDGId + p2.partonId.PDGId == -1) {
    emerging.partonId.PDGId = -24;
    emerging.colour = 0;
    emerging.antiColour = 0;
    isEWK = true;
    return emerging;
  }

  isEWK = false;
  return emerging;

}

vector<ClusteringConfigurationPtr> 
EWKClusterer::configurations (const vector<ClusteringParticleData>& in) {
  bool isEWK = false;
  vector<ClusteringConfigurationPtr> tmp;
  ClusteringParticleData emerging = emergingLine(in[0],in[1],isEWK);
  if (isEWK) {
    vector<ClusteringParticleData> out; out.push_back(emerging);
    tmp.push_back(new_ptr(ClusteringConfiguration(in,out,ClusteringInteractionType::EWK,this)));
  }
  return tmp;
}

void EWKClusterer::doKinematics (const tClusteringPtr& clustering) {
  Lorentz5Momentum emerging = clustering->children()[0]->momentum() + clustering->children()[1]->momentum();
  clustering->parents()[0]->momentum(emerging);

  Energy2 scale = clustering->scale();

  clustering->children()[0]->productionScale(scale);
  clustering->children()[1]->productionScale(scale);
  clustering->parents()[0]->splittingScale(scale);

  // the clustering has been chosen to be performed,
  // hoever, if we have already reached the electroweak
  // order, we don't consider it
  if (!performEWK()) clustering->veto();

}

ClusteringPtr EWKClusterer::doScale (const vector<tClusteringParticlePtr>& children,
					     const vector<ClusteringParticlePtr>& parents,
					     const tClusteringConfigurationPtr& config) {

  ClusteringPtr clustering = new_ptr(Clustering(children,parents,this,config));
  clustering->eventGenerator(generator());

  Energy2 scale = jetMeasure()->scale(children,parents,config);
  double z = jetMeasure()->z(children,parents,config);

  clustering->scale(scale);
  clustering->alphaScale(scale);
  clustering->momentumFraction(z);

  clustering->weight(1.);

  return clustering;

}
