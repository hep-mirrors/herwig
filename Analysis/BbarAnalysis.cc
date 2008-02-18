// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BbarAnalysis class.
//

#include "BbarAnalysis.h"

#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Histogram.h"


using namespace Herwig;
using namespace ThePEG;
using namespace std;

Histogram cos3_h(-1.,1.,100)     ,cos4_h(-1.,1.,100);
Histogram eta3_h(-3.,3.,100)     ,pt3_h(0.,100.,100);
Histogram eta4_h(-3.,3.,100)     ,pt4_h(0.,100.,100);
Histogram eta5_h(-3.,3.,100)     ,pt5_h(0.,100.,100);
Histogram eta34_h(-3.,3.,100)    ,y34_h(-4.,4.,80);
Histogram m34_h(0.,200.,40)     ,m345_h(0.,200.,40);

BbarAnalysis::~BbarAnalysis() {}


void BbarAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze( event, ieve, loop, state);
 
}

LorentzRotation BbarAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void BbarAnalysis::analyze(const tPVector & particles) {
  Lorentz5Momentum elec;
  Lorentz5Momentum pos;


  //find and store all FS electron / positron momenta
  for ( int i = 0; i <  particles.size(); ++i ){
    if( particles[i]->id() > 10 && particles[i]->id() < 16) elec = particles[i]->momentum();
    if( particles[i]->id() < -10 && particles[i]->id() > -16 ) pos = particles[i]->momentum();
  }

  Energy2 _mll2 = ( pos + elec ).m2();
  Energy m345 = 0. * GeV;
  double y34 = ( pos + elec ).rapidity();



  cos3_h.addWeighted(elec.cosTheta(),1.);
  cos4_h.addWeighted(pos.cosTheta(),1.);
  eta3_h.addWeighted(elec.eta(),1.);
  eta4_h.addWeighted(pos.eta(),1.);
  pt3_h.addWeighted(elec.perp(),1.);
  pt4_h.addWeighted(pos.perp(),1.);
  y34_h.addWeighted(y34,1.);
  m34_h.addWeighted(sqrt(_mll2)/GeV,1.);
  m345_h.addWeighted(sqrt(m345),1.);

  AnalysisHandler::analyze(particles);
}

void BbarAnalysis::analyze(tPPtr part) {
  //find electron (and positron) with highest pt
  //check to see this is the same as the Z pt
  //   if( part->id() == ParticleID::eminus ) {
  //  cerr<< "id is: "<<part->id()<<" pt = "<< part->momentum().perp() / GeV<<" \n ";
  // }
}

void BbarAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void BbarAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<BbarAnalysis> BbarAnalysis::initBbarAnalysis;
// Definition of the static class description member.

void BbarAnalysis::Init() {

  static ClassDocumentation<BbarAnalysis> documentation
    ("There is no documentation for the BbarAnalysis class");

}

void BbarAnalysis::dofinish() {
  AnalysisHandler::dofinish();

  //normalise all these histograms

  y34_h.normaliseToCrossSection();
  m34_h.normaliseToCrossSection();
  cos3_h.normaliseToCrossSection();
  cos4_h.normaliseToCrossSection();
  eta3_h.normaliseToCrossSection();
  eta4_h.normaliseToCrossSection();
  pt3_h.normaliseToCrossSection();
  pt4_h.normaliseToCrossSection();
  m345_h.normaliseToCrossSection();

  ofstream outfile("BbarAnalysis_hists.top");
  using namespace HistogramOptions;

  cos3_h.topdrawOutput(outfile,Frame,"RED","cos3 distribution: all wgts");
 
  cos4_h.topdrawOutput(outfile,Frame,"RED","cos4 distribution: all wgts");

  eta3_h.topdrawOutput(outfile,Frame,"RED","eta3 distribution: all wgts");
 
  eta4_h.topdrawOutput(outfile,Frame,"RED","eta4 distribution: all wgts");
 
  pt3_h.topdrawOutput(outfile,Frame,"RED","pt3 distribution: all wgts");
 
  pt4_h.topdrawOutput(outfile,Frame,"RED","pt4 distribution: all wgts");
 

  y34_h.topdrawOutput(outfile,Frame,"RED","y34 distribution: all wgts");
 
  m34_h.topdrawOutput(outfile,Frame,"RED","m34 distribution: all wgts");


  m345_h.topdrawOutput(outfile,Frame,"RED","m345 distribution: all wgts");

  outfile.close();
}
