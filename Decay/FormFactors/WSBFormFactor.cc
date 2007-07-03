// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WSBFormFactor class.
//

#include "WSBFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WSBFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;

WSBFormFactor::WSBFormFactor() 
{
  // modes handled by this and the parameters model
  // K to pi
  _F0.push_back(0.992);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(0.494*GeV);_mV0.push_back(0.892*GeV);
  _mS1.push_back(1.430*GeV);_mV1.push_back(1.273*GeV);
  _F0.push_back(0.992);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(0.498*GeV);_mV0.push_back(0.896*GeV);
  _mS1.push_back(1.430*GeV);_mV1.push_back(1.273*GeV);
  addFormFactor(321, 111,0,2,-3,-2);
  addFormFactor(311,-211,0,1,-3,-2);
  // D to K
  _F0.push_back(0.762);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  _F0.push_back(0.762);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  // D to pi
  _F0.push_back(0.692);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.692);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.692);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.692);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  // D to eta
  _F0.push_back(0.681);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(421,221,0,-2,4,2);
  _F0.push_back(0.681);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(411,221,0,-1,4,1);
  // D to eta'
  _F0.push_back(0.655);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(421,331,0,-2,4,2);
  _F0.push_back(0.655);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(411,331,0,-1,4,1);
  // D to K*
  _F0.push_back(0.000);_V.push_back(1.226);_A0.push_back(0.733);
  _A1.push_back(0.880);_A2.push_back(1.147);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  _F0.push_back(0.000);_V.push_back(1.226);_A0.push_back(0.733);
  _A1.push_back(0.880);_A2.push_back(1.147);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  // D to rho
  _F0.push_back(0.000);_V.push_back(1.225);_A0.push_back(0.669);
  _A1.push_back(0.775);_A2.push_back(0.923);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.000);_V.push_back(1.225);_A0.push_back(0.669);
  _A1.push_back(0.775);_A2.push_back(0.923);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.000);_V.push_back(1.225);_A0.push_back(0.669);
  _A1.push_back(0.775);_A2.push_back(0.923);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  _F0.push_back(0.000);_V.push_back(1.225);_A0.push_back(0.669);
  _A1.push_back(0.775);_A2.push_back(0.923);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  // D to omega
  _F0.push_back(0.000);_V.push_back(1.236);_A0.push_back(0.669);
  _A1.push_back(0.772);_A2.push_back(0.920);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(411,223,1,-1,4,1);
  _F0.push_back(0.000);_V.push_back(1.236);_A0.push_back(0.669);
  _A1.push_back(0.772);_A2.push_back(0.920);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(421,223,1,-2,4,2);
  // D_s to eta
  _F0.push_back(0.723);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  addFormFactor(431,221,0,-3,4,3);
  // D_s to eta'
  _F0.push_back(0.704);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  addFormFactor(431,331,0,-3,4,3);
  // D_s to K
  _F0.push_back(0.643);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(431,311,0,-3,4,1);
  _F0.push_back(0.643);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(431,321,0,-3,4,2);
  // D_s to K*
  _F0.push_back(0.000);_V.push_back(1.250);_A0.push_back(0.634);
  _A1.push_back(0.717);_A2.push_back(0.853);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(431,313,1,-3,4,1);
  _F0.push_back(0.000);_V.push_back(1.250);_A0.push_back(0.634);
  _A1.push_back(0.717);_A2.push_back(0.853);
  _mS0.push_back(1.87*GeV);_mV0.push_back(2.01*GeV);
  _mS1.push_back(2.47*GeV);_mV1.push_back(2.42*GeV);
  addFormFactor(431,323,1,-3,4,2);
  // D_s to phi
  _F0.push_back(0.000);_V.push_back(1.319);_A0.push_back(0.700);
  _A1.push_back(0.820);_A2.push_back(1.076);
  _mS0.push_back(1.97*GeV);_mV0.push_back(2.11*GeV);
  _mS1.push_back(2.60*GeV);_mV1.push_back(2.53*GeV);
  addFormFactor(431,333,1,-3,4,3);
  // B to D
  _F0.push_back(0.690);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(6.30*GeV);_mV0.push_back(6.34*GeV);
  _mS1.push_back(6.80*GeV);_mV1.push_back(6.73*GeV);
  _F0.push_back(0.690);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(6.30*GeV);_mV0.push_back(6.34*GeV);
  _mS1.push_back(6.80*GeV);_mV1.push_back(6.73*GeV);
  addFormFactor(521,-421,0,2,-5,-4);
  addFormFactor(511,-411,0,2,-5,-4);
  // B to K 
  _F0.push_back(0.379);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.38*GeV);_mV0.push_back(5.43*GeV);
  _mS1.push_back(5.89*GeV);_mV1.push_back(5.82*GeV);
  _F0.push_back(0.379);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.38*GeV);_mV0.push_back(5.43*GeV);
  _mS1.push_back(5.89*GeV);_mV1.push_back(5.82*GeV);
  addFormFactor(521,321,0,2,-5,-3);
  addFormFactor(511,311,0,1,-5,-3);
  // B to pi
  _F0.push_back(0.333);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  _F0.push_back(0.333);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521, 111,0,2,-5,-2);
  addFormFactor(511,-211,0,1,-5,-2);
  _F0.push_back(0.333);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  _F0.push_back(0.333);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521, 211,0,2,-5,-1);
  addFormFactor(511, 111,0,1,-5,-1);
  // B to eta
  _F0.push_back(0.307);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521,221,0,2,-5,-2);
  _F0.push_back(0.307);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(511,221,0,1,-5,-1);
  // B to eta'
  _F0.push_back(0.254);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521,331,0,2,-5,2);
  _F0.push_back(0.254);_V.push_back(0.000);_A0.push_back(0.000);
  _A1.push_back(0.000);_A2.push_back(0.000);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(511,331,0,1,-5,1);
  // B to D*
  _F0.push_back(0.000);_V.push_back(0.705);_A0.push_back(0.623);
  _A1.push_back(0.651);_A2.push_back(0.686);
  _mS0.push_back(6.30*GeV);_mV0.push_back(6.34*GeV);
  _mS1.push_back(6.80*GeV);_mV1.push_back(6.73*GeV);
  _F0.push_back(0.000);_V.push_back(0.705);_A0.push_back(0.623);
  _A1.push_back(0.651);_A2.push_back(0.686);
  _mS0.push_back(6.30*GeV);_mV0.push_back(6.34*GeV);
  _mS1.push_back(6.80*GeV);_mV1.push_back(6.73*GeV);
  addFormFactor(521,-423,1,2,-5,-4);
  addFormFactor(511,-413,1,1,-5,-4);
  // B to K*
  _F0.push_back(0.000);_V.push_back(0.369);_A0.push_back(0.321);
  _A1.push_back(0.328);_A2.push_back(0.331);
  _mS0.push_back(5.38*GeV);_mV0.push_back(5.43*GeV);
  _mS1.push_back(5.89*GeV);_mV1.push_back(5.82*GeV);
  _F0.push_back(0.000);_V.push_back(0.369);_A0.push_back(0.321);
  _A1.push_back(0.328);_A2.push_back(0.331);
  _mS0.push_back(5.38*GeV);_mV0.push_back(5.43*GeV);
  _mS1.push_back(5.89*GeV);_mV1.push_back(5.82*GeV);
  addFormFactor(521,323,1,2,-5,-3);
  addFormFactor(511,313,1,1,-5,-3);
  // B to rho 
  _F0.push_back(0.000);_V.push_back(0.329);_A0.push_back(0.281);
  _A1.push_back(0.283);_A2.push_back(0.283);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  _F0.push_back(0.000);_V.push_back(0.329);_A0.push_back(0.281);
  _A1.push_back(0.283);_A2.push_back(0.283);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521, 113,1,2,-5,-2);
  addFormFactor(511,-213,1,1,-5,-2);
  _F0.push_back(0.000);_V.push_back(0.329);_A0.push_back(0.281);
  _A1.push_back(0.283);_A2.push_back(0.283);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  _F0.push_back(0.000);_V.push_back(0.329);_A0.push_back(0.281);
  _A1.push_back(0.283);_A2.push_back(0.283);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521, 213,1,2,-5,-1);
  addFormFactor(511, 113,1,1,-5,-1);
  // B to omega
  _F0.push_back(0.000);_V.push_back(0.328);_A0.push_back(0.280);
  _A1.push_back(0.281);_A2.push_back(0.281);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(521,223,1,2,-5,-2);
  _F0.push_back(0.000);_V.push_back(0.328);_A0.push_back(0.280);
  _A1.push_back(0.281);_A2.push_back(0.281);
  _mS0.push_back(5.27*GeV);_mV0.push_back(5.32*GeV);
  _mS1.push_back(5.78*GeV);_mV1.push_back(5.71*GeV);
  addFormFactor(511,223,1,1,-5,-1);
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta=-0.194;
}

WSBFormFactor::~WSBFormFactor() {}

void WSBFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _F0 << _V << _A0 << _A1 << _A2 
     << ounit(_mS0,GeV) << ounit(_mS1,GeV) << ounit(_mV0,GeV) << ounit(_mV1,GeV) << _thetaeta;
}

void WSBFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _F0 >> _V >> _A0 >> _A1 >> _A2 
     >> iunit(_mS0,GeV) >> iunit(_mS1,GeV) >> iunit(_mV0,GeV) >> iunit(_mV1,GeV) >> _thetaeta;
}

ClassDescription<WSBFormFactor> WSBFormFactor::initWSBFormFactor;
// Definition of the static class description member.

void WSBFormFactor::Init() {

  static ClassDocumentation<WSBFormFactor> documentation
    ("The WSBFormFactor class is the implementation of the form-factors of "
     "Z.Phys.C29,637.");

  static ParVector<WSBFormFactor,double> interfaceF0
    ("F0",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_F0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceV
    ("V",
     "The form-factor V at zero q^2",
     &WSBFormFactor::_V,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA0
    ("A0",
     "The form-factor A0 at zero q^2",
     &WSBFormFactor::_A0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA1
    ("A1",
     "The form-factor A1 at zero q^2",
     &WSBFormFactor::_A1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA2
    ("A2",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_A2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceScalarMass
    ("ScalarMass",
     "The scalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS0,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoScalarMass
    ("PseudoScalarMass",
     "The pseudoscalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceVectorMass
    ("VectorMass",
     "The vector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV0,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoVectorMass
    ("PseudoVectorMass",
     "The pseudovector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static Parameter<WSBFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &WSBFormFactor::_thetaeta, -0.194, -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void WSBFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
					   int,int id1,
					   Energy, Energy,Complex & f0,
					   Complex & fp) const
{
  f0 = _F0[mode]/(1.-q2/_mS1[mode]/_mS1[mode]);
  fp = _F0[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect))
    {
      double fact;
      if(id1==ParticleID::eta)
	{
	  if(abs(outquark)==3){fact=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
	  else{fact=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
	}
      else if(id1==ParticleID::etaprime)
	{
	  if(abs(outquark)==3){fact=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
	  else{fact=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
	}
      else if(id1==ParticleID::pi0&&abs(outquark)==1){fact=-sqrt(0.5);}
      else{fact= sqrt(0.5);}
      f0*=fact;fp*=fact;
    }
}

void WSBFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
					   int, int id1, 
					   Energy, Energy,Complex & A0,
					   Complex & A1,Complex & A2,Complex & V) const
{
  A0 = -_A0[mode]/(1.-q2/_mS0[mode]/_mS0[mode]);
  A1 = -_A1[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  A2 = -_A2[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  V  = - _V[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3)
    {
      double fact(sqrt(0.5));
      if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
      A0*=fact;A1*=fact;A2*=fact;V*=fact;
    }
}

void WSBFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::WSBFormFactor " << fullName() << " \n";}
  output << "set " << fullName() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "set " << fullName() << ":F0 " 
		 << ix << "  " << _F0[ix] << endl;
	  output << "set " << fullName() << ":V  " 
		 << ix << "  " << _V[ix]  << endl;
	  output << "set " << fullName() << ":A0 " 
		 << ix << "  " << _A0[ix] << endl;
	  output << "set " << fullName() << ":A1 " 
		 << ix << "  " << _A1[ix] << endl;
	  output << "set " << fullName() << ":A2 " 
		 << ix << "  " << _A2[ix] << endl;
	  output << "set " << fullName() << ":ScalarMass " 
		 << ix << "  " << _mS0[ix]/GeV << endl;
	  output << "set " << fullName() << ":PseudoScalarMass " 
		 << ix << "  " << _mS1[ix]/GeV << endl;
	  output << "set " << fullName() << ":VectorMass " 
		 << ix << "  " << _mV0[ix]/GeV << endl;
	  output << "set " << fullName() << ":PseudoVectorMass " 
		 << ix << "  " << _mV1[ix]/GeV << endl;
	}
      else
	{
	  output << "insert " << fullName() << ":F0 " 
		 << ix << "  " << _F0[ix] << endl;
	  output << "insert " << fullName() << ":V  " 
		 << ix << "  " << _V[ix]  << endl;
	  output << "insert " << fullName() << ":A0 " 
		 << ix << "  " << _A0[ix] << endl;
	  output << "insert " << fullName() << ":A1 " 
		 << ix << "  " << _A1[ix] << endl;
	  output << "insert " << fullName() << ":A2 " 
		 << ix << "  " << _A2[ix] << endl;
	  output << "insert " << fullName() << ":ScalarMass " 
		 << ix << "  " << _mS0[ix]/GeV << endl;
	  output << "insert " << fullName() << ":PseudoScalarMass " 
		 << ix << "  " << _mS1[ix]/GeV << endl;
	  output << "insert " << fullName() << ":VectorMass " 
		 << ix << "  " << _mV0[ix]/GeV << endl;
	  output << "insert " << fullName() << ":PseudoVectorMass " 
		 << ix << "  " << _mV1[ix]/GeV << endl;
	}
    }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
