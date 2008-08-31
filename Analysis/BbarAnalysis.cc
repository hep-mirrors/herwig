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

Histogram cos3_h(-1., 1.,40)    , cos4_h(-1., 1.,40)  ;
Histogram HT34_h(0.,500.,25)    ;
Histogram eta3_h(-6.,6.,100)    , pt3_h(0.,200.,100)  ;
Histogram eta4_h(-6.,6.,100)    , pt4_h(0.,200.,100)  ;
Histogram eta34_h(-6.,6.,100)   , y34_h(-5.,5.,100)   , pt34_h(0.,500.,25);
Histogram m34wz_h(70.,100.,100) , m34h_h(114.95,115.05,100 );
Histogram eta5_h(-6.,6.,100)    , pt5_h(0.,200.,100)  ;
Histogram eta6_h(-6.,6.,100)    , pt6_h(0.,200.,100)  ;
Histogram eta56_h(-6.,6.,100)   , y56_h(-5.,5.,100)   , pt56_h(0.,150.,30);
Histogram m56_h(114.95,115.05,100)  ;
Histogram y3456_h(-4.,4.,80)    , m3456_h(0.,1000.,100);
Histogram HT_h(0.,500.,25)      ;

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
  Lorentz5Momentum pgluon   ;
  tPPtr l, lbar;
  tPPtr b, bbar;
  tPPtr gluon;

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
    if(particles[i]->id()== 5) b     = particles[i];
    if(particles[i]->id()==-5) bbar  = particles[i];
    if(particles[i]->id()==21) gluon = particles[i];
  }

  // Forbid analysis of events containing only a single lepton:
  if((l&&!lbar)||(lbar&&!l)) throw Exception() 
    << "BbarAnalysis::analyze\n"
    << "Cannot have just one lepton." << Exception::abortnow;

  if(l)     pl     = l->momentum()    ;
  if(lbar)  plbar  = lbar->momentum() ;
  if(b   )  pb     = b->momentum()    ;
  if(bbar)  pbbar  = bbar->momentum() ;
  if(gluon) pgluon = gluon->momentum();

  Lorentz5Momentum pV, pH, pVStar;
  pV     = pl+plbar;
  pH     = pb+pbbar;
  pVStar = pV+pH   ;
  Lorentz5Momentum pVrest(0.*GeV,0.*GeV,0.*GeV,pV.m(),pV.m());
  Lorentz5Momentum pl_Vrest, plbar_Vrest;
  pl_Vrest    = boostx(pl   , pV, pVrest);
  plbar_Vrest = boostx(plbar, pV, pVrest);

  // Lepton:
  cos3_h.addWeighted(pl_Vrest.cosTheta(),1.);
  // Anti-lepton:
  cos4_h.addWeighted(plbar_Vrest.cosTheta(),1.);
  // Scalar sum of lepton pts:
  HT34_h.addWeighted((pl.perp()+plbar.perp())/GeV,1.);
  // Lepton:
  eta3_h.addWeighted(pl.eta(),1.);
  pt3_h.addWeighted(pl.perp()/GeV,1.);
  // Anti-lepton:
  eta4_h.addWeighted(plbar.eta(),1.);
  pt4_h.addWeighted(plbar.perp()/GeV,1.);
  // Vector boson or Higgs boson depending on process:
  eta34_h.addWeighted(pV.eta(),1.);
  y34_h.addWeighted(pV.rapidity(),1.);
  pt34_h.addWeighted(pV.perp()/GeV,1.);
  m34wz_h.addWeighted(sqrt(pV.m2())/GeV,1.);
  m34h_h.addWeighted(sqrt(pV.m2())/GeV,1.);
  // b-quark:
  eta5_h.addWeighted(pb.eta(),1.);
  pt5_h.addWeighted(pb.perp()/GeV,1.);
  // Anti-b-quark:
  eta6_h.addWeighted(pbbar.eta(),1.);
  pt6_h.addWeighted(pbbar.perp()/GeV,1.);
  // Higgs bosons only:
  y56_h.addWeighted(pH.rapidity(),1.);
  eta56_h.addWeighted(pH.eta(),1.);
  pt56_h.addWeighted(pH.perp()/GeV,1.);
  m56_h.addWeighted(sqrt(pH.m2())/GeV,1.);
  // Everything added up:
  y3456_h.addWeighted(pVStar.rapidity(),1.);
  m3456_h.addWeighted(sqrt(pVStar.m2())/GeV,1.);


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
  // Anti-lepton:
  cos4_h.normaliseToCrossSection();
  cos4_h.prefactor(cos4_h.prefactor()*1.e6);
  // HT34:
  HT34_h.normaliseToCrossSection();
  HT34_h.prefactor(HT34_h.prefactor()*1.e6);
  // Lepton:
  eta3_h.normaliseToCrossSection();
  eta3_h.prefactor(eta3_h.prefactor()*1.e6);
  pt3_h.normaliseToCrossSection();
  pt3_h.prefactor(pt3_h.prefactor()*1.e6);
  // Anti-lepton:
  eta4_h.normaliseToCrossSection();
  eta4_h.prefactor(eta4_h.prefactor()*1.e6);
  pt4_h.normaliseToCrossSection();
  pt4_h.prefactor(pt4_h.prefactor()*1.e6);
  // Vector boson or Higgs boson depending on process:
  eta34_h.normaliseToCrossSection();
  eta34_h.prefactor(eta34_h.prefactor()*1.e6);
  y34_h.normaliseToCrossSection();
  y34_h.prefactor(y34_h.prefactor()*1.e6);
  pt34_h.normaliseToCrossSection();
  pt34_h.prefactor(pt34_h.prefactor()*1.e6);
  m34wz_h.normaliseToCrossSection();
  m34wz_h.prefactor(m34wz_h.prefactor()*1.e6);
  m34h_h.normaliseToCrossSection();
  m34h_h.prefactor(m34h_h.prefactor()*1.e6);
  // b-quark:
  eta5_h.normaliseToCrossSection();
  eta5_h.prefactor(eta5_h.prefactor()*1.e6);
  pt5_h.normaliseToCrossSection();
  pt5_h.prefactor(pt5_h.prefactor()*1.e6);
  // Anti-b-quark:
  eta6_h.normaliseToCrossSection();
  eta6_h.prefactor(eta6_h.prefactor()*1.e6);
  pt6_h.normaliseToCrossSection();
  pt6_h.prefactor(pt6_h.prefactor()*1.e6);
  // Higgs bosons only:
  y56_h.normaliseToCrossSection();
  y56_h.prefactor(y56_h.prefactor()*1.e6);
  eta56_h.normaliseToCrossSection();
  eta56_h.prefactor(eta56_h.prefactor()*1.e6);
  pt56_h.normaliseToCrossSection();
  pt56_h.prefactor(pt56_h.prefactor()*1.e6);
  m56_h.normaliseToCrossSection();
  m56_h.prefactor(m56_h.prefactor()*1.e6);
  // Everything added up:
  y3456_h.normaliseToCrossSection();
  y3456_h.prefactor(y3456_h.prefactor()*1.e6);
  m3456_h.normaliseToCrossSection();
  m3456_h.prefactor(m3456_h.prefactor()*1.e6);

  file.open(fname.c_str());
  using namespace HistogramOptions;

  // Lepton:
  cos3_h.topdrawOutput(file,Frame,"RED","cos3 distribution: all wgts");
  // Anti-lepton:
  cos4_h.topdrawOutput(file,Frame,"RED","cos4 distribution: all wgts");
  // HT34:
  HT34_h.topdrawOutput(file,Frame,"RED","HT34 distribution: all wgts");
  // Lepton:
  eta3_h.topdrawOutput(file,Frame,"RED","eta3 distribution: all wgts");
  pt3_h.topdrawOutput(file,Frame,"RED","pt3 distribution: all wgts");
  // Anti-lepton:
  eta4_h.topdrawOutput(file,Frame,"RED","eta4 distribution: all wgts");
  pt4_h.topdrawOutput(file,Frame,"RED","pt4 distribution: all wgts");
  // Vector boson or Higgs boson depending on process:
  eta34_h.topdrawOutput(file,Frame,"RED","eta34 distribution: all wgts");
  y34_h.topdrawOutput(file,Frame,"RED","y34 distribution: all wgts");
  pt34_h.topdrawOutput(file,Frame,"RED","pt34 low  mass: all wgts");
  m34wz_h.topdrawOutput(file,Frame,"RED","m34 low  mass: all wgts");
  m34h_h.topdrawOutput(file,Frame,"RED","m34 high mass: all wgts");
  // b-quark:
  eta5_h.topdrawOutput(file,Frame,"RED","eta5 distribution: all wgts");
  pt5_h.topdrawOutput(file,Frame,"RED","pt5 distribution: all wgts");
  // Anti-b-quark:
  eta6_h.topdrawOutput(file,Frame,"RED","eta6 distribution: all wgts");
  pt6_h.topdrawOutput(file,Frame,"RED","pt6 distribution: all wgts");
  // Higgs bosons only:
  y56_h.topdrawOutput(file,Frame,"RED","y56 distribution: all wgts");
  eta56_h.topdrawOutput(file,Frame,"RED","eta56 distribution: all wgts");
  pt56_h.topdrawOutput(file,Frame,"RED","eta56 distribution: all wgts");
  m56_h.topdrawOutput(file,Frame,"RED","m56 distribution: all wgts");
  // Everything added up:
  y3456_h.topdrawOutput(file,Frame,"RED","y3456 distribution: all wgts");
  m3456_h.topdrawOutput(file,Frame,"RED","m3456 distribution: all wgts");

  file.close();
}


Lorentz5Momentum BbarAnalysis::boostx(Lorentz5Momentum p_in, 
				      Lorentz5Momentum pt, 
				      Lorentz5Momentum ptt){
  // Boost input vector p_in to output vector p_out using the same
  // transformation as required to boost massive vector pt to ptt
  Lorentz5Momentum p_tmp,p_out;
  Vector3<double> beta;
  Energy  mass, bdotp;
  double  gam;

  if (pt.m2()<0.*GeV2) {
    throw Exception() << "pt.m2()<0. in boostx, pt.m2() = "
		      << pt.m2()/GeV2 
		      << Exception::runerror;
  }
  mass = sqrt(pt.m2());

  // boost to the rest frame of pt
  gam = pt.e()/mass;

  beta  = -pt.vect()*(1./pt.e());
  bdotp =  beta*p_in.vect();

  Lorentz5Momentum tmp_vec;
  tmp_vec.setVect(beta*gam*((p_in.e()+gam*(p_in.e()+bdotp))/(1.+gam)));
  tmp_vec.setE(0.*GeV);
  tmp_vec.setMass(0.*GeV);
  p_tmp = p_in + tmp_vec;
  p_tmp.setE(gam*(p_in.e()+bdotp));

  // boost from rest frame of pt to frame in which pt is identical
  // with ptt, thus completing the transformation
  gam   =  ptt.e()/mass;
  beta  =  ptt/ptt.e();
  bdotp =  beta*p_tmp.vect();

  p_out.setVect(p_tmp.vect()+gam*beta*((p_out.e()+p_tmp.e())/(1.+gam)));
  p_out.setE(gam*(p_tmp.e()+bdotp));
  p_out.rescaleMass();

  return p_out;
}
