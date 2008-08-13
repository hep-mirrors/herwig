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

Histogram cos3_h(-1.,1.,100),eta3_h(-6.,6.,120)   ,pt3_h(0.,200.,100);
Histogram cos4_h(-1.,1.,100),eta4_h(-6.,6.,120)   ,pt4_h(0.,200.,100);
Histogram cos5_h(-1.,1.,100),eta5_h(-6.,6.,120)   ,pt5_h(0.,200.,100);
Histogram cos6_h(-1.,1.,100),eta6_h(-6.,6.,120)   ,pt6_h(0.,200.,100);
Histogram y34_h(-6.,6.,120)  ,eta34_h(-6.,6.,120)  ,m34wz_h(70.,100.,100 );
Histogram                                          m34h_h(114.,116.,100 );
Histogram y56_h(-6.,6.,120)  ,eta56_h(-6.,6.,120)  ,m56_h(114.,116.,100 );
Histogram y3456_h(-6.,6.,120),eta3456_h(-6.,6.,120),m3456_h(100.,600.,5 );

BbarAnalysis::~BbarAnalysis() {}


void BbarAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze( event, ieve, loop, state);
 
}

LorentzRotation BbarAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void BbarAnalysis::analyze(const tPVector & particles) {
  Lorentz5Momentum pl, plbar;
  Lorentz5Momentum pb, pbbar;
  tPPtr l, lbar;
  tPPtr b, bbar;

  // Find and store all lepton and antilepton pairs
  // and likewise any b-bbar pairs.
  for (unsigned int i = 0; i <  particles.size(); ++i ){
    // Avoid anything that might've come from a remnant decay.
    if(particles[i]->parents().size()>0&&particles[i]->parents()[0]->id()==82)
      continue;
    if(particles[i]->id()> 10&&particles[i]->id()< 17) 
      l    = particles[i];
    if(particles[i]->id()<-10&&particles[i]->id()>-17) 
      lbar = particles[i];
    if(particles[i]->id()== 5) b    = particles[i];
    if(particles[i]->id()==-5) bbar = particles[i];
  }

  // Forbid analysis of events containing only a single lepton:
  if((l&&!lbar)||(lbar&&!l)) throw Exception() 
    << "BbarAnalysis::analyze\n"
    << "Cannot have just one lepton." << Exception::abortnow;

  if(l)    pl    = l->momentum()   ;
  if(lbar) plbar = lbar->momentum();
  if(b   ) pb    = b->momentum()   ;
  if(bbar) pbbar = bbar->momentum();

  // Lepton:
  cos3_h.addWeighted(pl.cosTheta(),1.);
  eta3_h.addWeighted(pl.eta(),1.);
  pt3_h.addWeighted(pl.perp()/GeV,1.);
  // Anti-lepton:
  cos4_h.addWeighted(plbar.cosTheta(),1.);
  eta4_h.addWeighted(plbar.eta(),1.);
  pt4_h.addWeighted(plbar.perp()/GeV,1.);
  // b-quark:
  cos5_h.addWeighted(pb.cosTheta(),1.);
  eta5_h.addWeighted(pb.eta(),1.);
  pt5_h.addWeighted(pb.perp()/GeV,1.);
  // Anti-b-quark:
  cos6_h.addWeighted(pbbar.cosTheta(),1.);
  eta6_h.addWeighted(pbbar.eta(),1.);
  pt6_h.addWeighted(pbbar.perp()/GeV,1.);
  // Vector boson or Higgs boson depending on process:
  y34_h.addWeighted((pl+plbar).rapidity(),1.);
  eta34_h.addWeighted((pl+plbar).eta(),1.);
  m34wz_h.addWeighted(sqrt((pl+plbar).m2())/GeV,1.);
  m34h_h.addWeighted(sqrt((pl+plbar).m2())/GeV,1.);
  // Higgs bosons only:
  y56_h.addWeighted((pb+pbbar).rapidity(),1.);
  eta56_h.addWeighted((pb+pbbar).eta(),1.);
  m56_h.addWeighted(sqrt((pb+pbbar).m2())/GeV,1.);
  // Everything added up:
  y3456_h.addWeighted((pl+plbar+pb+pbbar).rapidity(),1.);
  eta3456_h.addWeighted((pl+plbar+pb+pbbar).eta(),1.);
  m3456_h.addWeighted(sqrt((pl+plbar+pb+pbbar).m2())/GeV,1.);


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
  ofstream file;
  string fname = generator()->filename()+string("_")+name()+string(".top");

  // Normalise the histograms
  // Lepton:
  cos3_h.normaliseToCrossSection();
  cos3_h.prefactor(cos3_h.prefactor()*1.e6);
  eta3_h.normaliseToCrossSection();
  eta3_h.prefactor(eta3_h.prefactor()*1.e6);
  pt3_h.normaliseToCrossSection();
  pt3_h.prefactor(pt3_h.prefactor()*1.e6);
  // Anti-lepton:
  cos4_h.normaliseToCrossSection();
  cos4_h.prefactor(cos4_h.prefactor()*1.e6);
  eta4_h.normaliseToCrossSection();
  eta4_h.prefactor(eta4_h.prefactor()*1.e6);
  pt4_h.normaliseToCrossSection();
  pt4_h.prefactor(pt4_h.prefactor()*1.e6);
  // b-quark:
  cos5_h.normaliseToCrossSection();
  cos5_h.prefactor(cos5_h.prefactor()*1.e6);
  eta5_h.normaliseToCrossSection();
  eta5_h.prefactor(eta5_h.prefactor()*1.e6);
  pt5_h.normaliseToCrossSection();
  pt5_h.prefactor(pt5_h.prefactor()*1.e6);
  // Anti-b-quark:
  cos6_h.normaliseToCrossSection();
  cos6_h.prefactor(cos6_h.prefactor()*1.e6);
  eta6_h.normaliseToCrossSection();
  eta6_h.prefactor(eta6_h.prefactor()*1.e6);
  pt6_h.normaliseToCrossSection();
  pt6_h.prefactor(pt6_h.prefactor()*1.e6);
  // Vector boson or Higgs boson depending on process:
  y34_h.normaliseToCrossSection();
  y34_h.prefactor(y34_h.prefactor()*1.e6);
  eta34_h.normaliseToCrossSection();
  eta34_h.prefactor(eta34_h.prefactor()*1.e6);
  m34wz_h.normaliseToCrossSection();
  m34wz_h.prefactor(m34wz_h.prefactor()*1.e6);
  m34h_h.normaliseToCrossSection();
  m34h_h.prefactor(m34h_h.prefactor()*1.e6);
  // Higgs bosons only:
  y56_h.normaliseToCrossSection();
  y56_h.prefactor(y56_h.prefactor()*1.e6);
  eta56_h.normaliseToCrossSection();
  eta56_h.prefactor(eta56_h.prefactor()*1.e6);
  m56_h.normaliseToCrossSection();
  m56_h.prefactor(m56_h.prefactor()*1.e6);
  // Everything added up:
  y3456_h.normaliseToCrossSection();
  y3456_h.prefactor(y3456_h.prefactor()*1.e6);
  eta3456_h.normaliseToCrossSection();
  eta3456_h.prefactor(eta3456_h.prefactor()*1.e6);
  m3456_h.normaliseToCrossSection();
  m3456_h.prefactor(m3456_h.prefactor()*1.e6);

  file.open(fname.c_str());
  using namespace HistogramOptions;

  // Lepton:
  cos3_h.topdrawOutput(file,Frame,"RED","cos3 distribution: all wgts");
  eta3_h.topdrawOutput(file,Frame,"RED","eta3 distribution: all wgts");
  pt3_h.topdrawOutput(file,Frame,"RED","pt3 distribution: all wgts");
  // Anti-lepton:
  cos4_h.topdrawOutput(file,Frame,"RED","cos4 distribution: all wgts");
  eta4_h.topdrawOutput(file,Frame,"RED","eta4 distribution: all wgts");
  pt4_h.topdrawOutput(file,Frame,"RED","pt4 distribution: all wgts");
  // b-quark:
  cos5_h.topdrawOutput(file,Frame,"RED","cos5 distribution: all wgts");
  eta5_h.topdrawOutput(file,Frame,"RED","eta5 distribution: all wgts");
  pt5_h.topdrawOutput(file,Frame,"RED","pt5 distribution: all wgts");
  // Anti-b-quark:
  cos6_h.topdrawOutput(file,Frame,"RED","cos6 distribution: all wgts");
  eta6_h.topdrawOutput(file,Frame,"RED","eta6 distribution: all wgts");
  pt6_h.topdrawOutput(file,Frame,"RED","pt6 distribution: all wgts");
  // Vector boson or Higgs boson depending on process:
  y34_h.topdrawOutput(file,Frame,"RED","y34 distribution: all wgts");
  eta34_h.topdrawOutput(file,Frame,"RED","eta34 distribution: all wgts");
  m34wz_h.topdrawOutput(file,Frame,"RED","m34 low  mass: all wgts");
  m34h_h.topdrawOutput(file,Frame,"RED","m34 high mass: all wgts");
  // Higgs bosons only:
  y56_h.topdrawOutput(file,Frame,"RED","y56 distribution: all wgts");
  eta56_h.topdrawOutput(file,Frame,"RED","eta56 distribution: all wgts");
  m56_h.topdrawOutput(file,Frame,"RED","m56 distribution: all wgts");
  // Everything added up:
  y3456_h.topdrawOutput(file,Frame,"RED","y3456 distribution: all wgts");
  eta3456_h.topdrawOutput(file,Frame,"RED","eta3456 distribution: all wgts");
  m3456_h.topdrawOutput(file,Frame,"RED","m3456 distribution: all wgts");

  file.close();
}
