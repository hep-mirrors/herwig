// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SingleParticleAnalysis class.
//

#include "SingleParticleAnalysis.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SingleParticleAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SingleParticleAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
}

LorentzRotation SingleParticleAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void SingleParticleAnalysis::analyze(const tPVector & particles) {
  int Ncharged = 0; 
  for(unsigned int ix=0;ix<particles.size();++ix)
    {
      const Lorentz5Momentum p = particles[ix]->momentum();
      if(!particles[ix]->data().charged()) continue;
      Ncharged++;
      *_yT += abs(_shapes->yT(p));
      *_yS += abs(_shapes->yS(p));
      *_ptinT += _shapes->ptInT(p)/GeV;
      *_ptoutT += _shapes->ptOutT(p)/GeV;
      *_ptinS += _shapes->ptInS(p)/GeV;
      *_ptoutS += _shapes->ptOutS(p)/GeV;
    }
  // make sure the right bin is booked by subtracting a small amount. 
  *_nch += Ncharged-0.00001;
}

void SingleParticleAnalysis::analyze(tPPtr) {}

void SingleParticleAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void SingleParticleAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<SingleParticleAnalysis> SingleParticleAnalysis::initSingleParticleAnalysis;
// Definition of the static class description member.

void SingleParticleAnalysis::Init() {

  static ClassDocumentation<SingleParticleAnalysis> documentation
    ("There is no documentation for the SingleParticleAnalysis class");

  static Reference<SingleParticleAnalysis,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &SingleParticleAnalysis::_shapes, false, false, true, false, false);
}

