// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QCDClusterer class.
//

#include "QCDClusterer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QCDClusterer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/Shower/CKKW/Clustering/ToIncomingCMS.h"

using namespace Herwig;

QCDClusterer::~QCDClusterer() {}

void QCDClusterer::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _useColour;
}

void QCDClusterer::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _useColour;
}

ClusteringParticleData QCDClusterer::doEmergingLine
(ClusteringParticleData p1, ClusteringParticleData p2, bool& isQCD) const {

  ClusteringParticleData emerging;

  isQCD = false;

  bool colourOk = false;

  // qf / gf

  if (abs(p1.partonId.PDGId) < 7 && p2.partonId.PDGId == 21 &&
      p1.partonId.state == ClusteringParticleState::final &&
      p2.partonId.state == ClusteringParticleState::final) {

    if (_useColour) {
      if (p1.partonId.PDGId < 0 && p1.antiColour == p2.colour) {
	emerging.antiColour = p2.antiColour;
	colourOk = true;
      }
      if (p1.partonId.PDGId > 0 && p1.colour == p2.antiColour) {
	emerging.colour = p2.colour;
	colourOk = true;
      }

      if (!colourOk) { isQCD = false; return emerging; }
    }

    emerging.partonId.PDGId = p1.partonId.PDGId;
    emerging.partonId.state = ClusteringParticleState::final;

    isQCD = true;
    return emerging;

  }


  // gf / gf

  if (p1.partonId.PDGId == 21 && p2.partonId.PDGId == 21 &&
      p2.partonId.state == ClusteringParticleState::final &&
      p1.partonId.state == ClusteringParticleState::final) {

    if (_useColour) {

      if (p1.colour == p2.antiColour && p1.antiColour != p2.colour) {
	emerging.antiColour = p1.antiColour;
	emerging.colour = p2.colour;
	colourOk = true;
      }

      if (p2.colour == p1.antiColour && p2.antiColour != p1.colour) {
	emerging.colour = p1.colour;
	emerging.antiColour = p2.antiColour;
	colourOk = true;
      }

      if (p1.colour == p1.antiColour || p2.colour == p2.antiColour)
	colourOk = false;

      if (!colourOk) { isQCD = false; return emerging; }

    }

    emerging.partonId.PDGId = 21;
    emerging.partonId.state = ClusteringParticleState::final;
    isQCD = true;
    return emerging;

  }


  // qf / qbarf

  if (p1.partonId.PDGId + p2.partonId.PDGId == 0 && abs(p1.partonId.PDGId) < 7 &&
      p2.partonId.state == ClusteringParticleState::final &&
      p1.partonId.state == ClusteringParticleState::final) {
    
    if (_useColour) {
      
      if (p1.partonId.PDGId < 0 && p1.antiColour != p2.colour) {
	emerging.antiColour = p1.antiColour;
	emerging.colour = p2.colour;
	colourOk = true;
      }

      if (p2.partonId.PDGId < 0 && p2.antiColour != p1.colour) {
	emerging.antiColour = p2.antiColour;
	emerging.colour = p1.colour;
	colourOk = true;
      }

      if (!colourOk) { isQCD = false; return emerging; }

    }

    emerging.partonId.PDGId = 21;
    emerging.partonId.state = ClusteringParticleState::final;
    isQCD = true;
    return emerging;

  }


  // qi / gf

  if (abs(p1.partonId.PDGId) < 7 && p2.partonId.PDGId == 21 &&
      p1.partonId.state == ClusteringParticleState::initial &&
      p2.partonId.state == ClusteringParticleState::final) {
    
    if (_useColour) {

      if (p1.partonId.PDGId > 0 && p1.colour == p2.colour) {
	emerging.colour = p2.antiColour;
	colourOk = true;
      }

      if (p1.partonId.PDGId < 0 && p1.antiColour == p1.antiColour) {
	emerging.antiColour = p2.colour;
	colourOk = true;
      }

      if (!colourOk) { isQCD = false; return emerging; }

    }

    emerging.partonId.PDGId = p1.partonId.PDGId;
    emerging.partonId.state = ClusteringParticleState::initial;
    isQCD = true;
    return emerging;

  }

  // gi / gf

  if (p1.partonId.PDGId == 21 && p2.partonId.PDGId == 21 &&
      p1.partonId.state == ClusteringParticleState::initial &&
      p2.partonId.state == ClusteringParticleState::final) {

    if (_useColour) {

      if (p1.colour == p2.colour && p1.antiColour != p2.antiColour) {
	emerging.colour = p2.antiColour;
	emerging.antiColour = p1.antiColour;
	colourOk = true;
      }

      if (p1.antiColour == p2.antiColour && p1.colour != p2.colour) {
	emerging.antiColour = p2.colour;
	emerging.colour = p1.colour;
	colourOk = true;
      }

      if (p1.colour == p1.antiColour || p2.colour == p2.antiColour)
	colourOk = false;

      if (!colourOk) { isQCD = false; return emerging; }

    }

    emerging.partonId.PDGId = 21;
    emerging.partonId.state = ClusteringParticleState::initial;
    isQCD = true;
  
  }


  // gi / qf

  if (p1.partonId.PDGId == 21 && abs(p2.partonId.PDGId) < 7 &&
      p1.partonId.state == ClusteringParticleState::initial &&
      p2.partonId.state == ClusteringParticleState::final) {

    if (_useColour) {

      if (p2.partonId.PDGId > 0 && p1.colour == p2.colour) {
	emerging.antiColour = p1.antiColour;
	colourOk = true;
      }

      if (p2.partonId.PDGId < 0 && p1.antiColour == p2.antiColour) {
	emerging.colour = p1.colour;
	colourOk = true;
      }

      if (!colourOk) { isQCD = false; return emerging; }
    }

    emerging.partonId.PDGId = p2.partonId.PDGId;
    emerging.partonId.state = ClusteringParticleState::initial;
    isQCD = true;

  }

  // qi / qf

  if (p1.partonId.PDGId == p2.partonId.PDGId && abs(p1.partonId.PDGId) < 7 &&
      p1.partonId.state == ClusteringParticleState::initial &&
      p2.partonId.state == ClusteringParticleState::final) {

    if (_useColour) {

      if (p1.partonId.PDGId > 0 && p1.colour != p2.colour) {
	emerging.colour = p1.colour;
	emerging.antiColour = p2.colour;
	colourOk = true;
      }

      if (p1.partonId.PDGId < 0 && p1.antiColour != p2.antiColour) {
	emerging.antiColour = p1.antiColour;
	emerging.colour = p2.antiColour;
	colourOk = true;
      }

      if (!colourOk) { isQCD = false; return emerging; }
    }

    emerging.partonId.PDGId = 21;
    emerging.partonId.state = ClusteringParticleState::initial;
    isQCD = true;

  }

  return emerging;

}

bool QCDClusterer::colourConnected (ClusteringParticleData p1, ClusteringParticleData p2) {

  if (!_useColour) {
    // any QCD particle is colour connected
    return ((p1.partonId.PDGId == 21 || abs(p1.partonId.PDGId) <7)
      && (p2.partonId.PDGId == 21 || abs(p2.partonId.PDGId) <7));
  }

  if (p1.colour == p1.antiColour || p2.colour == p2.antiColour) return false;

  if (p1.partonId.state == p2.partonId.state) {
    return ((p1.colour == p2.antiColour && p1.colour !=0)
	    || (p1.antiColour == p2.colour && p1.antiColour != 0));
  }

  if (p1.partonId.state != p2.partonId.state) {
    return ((p1.colour == p2.colour && p1.colour != 0)
	    || (p1.antiColour == p2.antiColour && p1.antiColour != 0));
  }

  return false;

}

vector<ClusteringConfigurationPtr> 
QCDClusterer::configurations (const vector<ClusteringParticleData>& in) {
  bool isQCD = false;
  vector<ClusteringConfigurationPtr> tmp;
  ClusteringParticleData emerging = emergingLine(in[0],in[1],isQCD);
  if (isQCD) {
    vector<ClusteringParticleData> out; out.push_back(emerging);
    tmp.push_back(new_ptr(ClusteringConfiguration(in,out,ClusteringInteractionType::QCD,this)));
  }
  return tmp;
}

void QCDClusterer::doKinematics (const tClusteringPtr& clustering) {
  Lorentz5Momentum emerging;
  bool initial = false;
  if (clustering->children()[0]->pData().partonId.state == ClusteringParticleState::final &&
      clustering->children()[1]->pData().partonId.state == ClusteringParticleState::final)
    emerging = clustering->children()[0]->momentum() + clustering->children()[1]->momentum();
  if (clustering->children()[0]->pData().partonId.state == ClusteringParticleState::initial &&
      clustering->children()[1]->pData().partonId.state == ClusteringParticleState::final) {
    emerging = clustering->children()[0]->momentum() - clustering->children()[1]->momentum();
    initial = true;
  }
  if (clustering->children()[0]->pData().partonId.state == ClusteringParticleState::final &&
      clustering->children()[1]->pData().partonId.state == ClusteringParticleState::initial) {
    emerging = - clustering->children()[0]->momentum() + clustering->children()[1]->momentum();
    initial = true;
  }
  clustering->parents()[0]->momentum(emerging);

  Energy2 scale = clustering->scale();

  clustering->children()[0]->productionScale(scale);
  clustering->children()[1]->productionScale(scale);
  clustering->parents()[0]->splittingScale(scale);

}

ClusteringPtr QCDClusterer::doScale (const vector<tClusteringParticlePtr>& children,
					     const vector<ClusteringParticlePtr>& parents,
					     const tClusteringConfigurationPtr& config) {

  // by default, we set the alpha_s scale to the clustering scale,
  // as well as the adjacent particle's splitting/production scales
  // the weight is chosen to be alpha_s at the clustering scale
  // the emerging initial parton's longitudinal momentum fraction is 
  // set to z x , where x is the clustered one's fraction.
  // if initial state partons are involved, a ToIncomingCMS object
  // is set as PostClustering

  ClusteringPtr clustering = new_ptr(Clustering(children,parents,this,config));
  clustering->eventGenerator(generator());

  Energy2 scale = jetMeasure()->scale(children,parents,config);
  double z = jetMeasure()->z(children,parents,config);

  clustering->scale(scale);
  clustering->alphaScale(scale);
  clustering->momentumFraction(z);

  if(parents[0]->pData().partonId.state == ClusteringParticleState::initial) {

    parents[0]->x(z * parents[0]->x());
    clustering->postClustering(new_ptr(ToIncomingCMS()));

  }

  clustering->weight(1.);

  return clustering;

}

ClassDescription<QCDClusterer> QCDClusterer::initQCDClusterer;
// Definition of the static class description member.

void QCDClusterer::Init() {

  static ClassDocumentation<QCDClusterer> documentation
    ("QCDClusterer provides mapping of quantum numbers "
     "for simple 2->1 clusterings according to QCD vertices.");

  static Switch<QCDClusterer,bool> interfaceUseColour
    ("UseColour",
     "Wether or not to use colour information on clusterings",
     &QCDClusterer::_useColour, true, false, false);
  static SwitchOption interfaceUseColourUseColourOn
    (interfaceUseColour,
     "Yes",
     "Do use colour information on clustering.",
     true);
  static SwitchOption interfaceUseColourUseColourOff
    (interfaceUseColour,
     "No",
     "Do not use colour information on clustering.",
     false);


}

