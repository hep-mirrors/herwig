// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VV_ME_Analysis class.
//

#include "VV_ME_Analysis.h"

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
Histogram m34wz_h(70.,100.,100) ;
Histogram eta5_h(-6.,6.,100)    , pt5_h(0.,200.,100)  ;
Histogram eta6_h(-6.,6.,100)    , pt6_h(0.,200.,100)  ;
Histogram eta56_h(-6.,6.,100)   , y56_h(-5.,5.,100)   , pt56_h(0.,150.,30);
Histogram m56_h(75.0,100.0,100) ;
Histogram HT3456_h(0.,500.,25)  ;
Histogram y3456_h(-4.,4.,80)    , m3456_h(0.,1000.,100);
Histogram cth1_h(-1.,1.,40)     , acth2_h(0.,Constants::pi,40);

VV_ME_Analysis::~VV_ME_Analysis() {}

void VV_ME_Analysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  // Store the initial state particles.
  inbound=event->incoming();
  // Store all of the final state particles as leptons.
  leptons=event->getFinalState();
  // Extract the emitted parton
  tPVector::iterator pit;
  int nemitted(0);
  for(pit=leptons.begin();pit!=leptons.end();++pit) {
    if(!emitted) {
      if(((*pit)->id()==21)||(abs((*pit)->id())>=1&&abs((*pit)->id()<=6))) { 
	emitted = *pit;
	nemitted++;
      }
    }
  }
  if(nemitted>1)
	throw Exception() << "VV_ME_Analysis::analyze"  
			  << "\nMore than one QCD charge in the final state"
			  << "\nCan't tell what the NLO emitted parton is."
			  << Exception::warning;
  // Now delete everything else in leptons which is not a lepton!
  for(pit=leptons.begin();pit!=leptons.end();++pit)
    if(!(abs((*pit)->id())>=11&&abs((*pit)->id())<=16)) {
      leptons.erase(pit);
      --pit;
    }
  // Make sure that leptons has only 4 entries (2 bosons * 2-2-body decays).
  assert(leptons.size()==4);  
  // Make sure all final-state leptons have parents.
  for (unsigned int i=0; i<leptons.size(); i++)
    assert(leptons[i]->parents().size()>0);
  // Go get the parents of the leptons W+ W- / W+ Z / W- Z / Z Z
  tPVector bosons;
  tPVector::iterator qit;
  bosons.push_back(leptons[0]->parents()[0]);
  for(pit=leptons.begin();pit!=leptons.end();++pit) {
    bool already_got_parent(0);
    for(qit=bosons.begin();qit!=bosons.end();++qit)
      if((*pit)->parents()[0]==(*qit)) already_got_parent = 1;
    if(!already_got_parent) bosons.push_back((*pit)->parents()[0]);
  }
  // Attempt to order bosons according mcfm notation.
  if(abs(bosons[0]->id())==abs(bosons[1]->id())) {
    // Either W W or Z Z here. Swap so that the W+ is the first entry.
    if(bosons[0]->id()==-24) swap(bosons[0],bosons[1]);
    if(bosons[0]->id()== 23 &&
       (abs(bosons[1]->children()[0]->id())==13 ||
	abs(bosons[1]->children()[0]->id())==14 )) 
      swap(bosons[0],bosons[1]);
  } else {
    // Either W+ Z or W- Z here. Swap so that the W is first entry.
    if(bosons[0]->id()== 23) swap(bosons[0],bosons[1]);
  }
  // Boson children are ordered fermion, anti-fermion in the mcfm convention.
  if(bosons[0]->children()[0]->id()>0) {
    leptons[0]=bosons[0]->children()[0];
    leptons[1]=bosons[0]->children()[1];
  } else {
    leptons[0]=bosons[0]->children()[1];
    leptons[1]=bosons[0]->children()[0];
  }
  if(bosons[1]->children()[0]->id()>0) {
    leptons[2]=bosons[1]->children()[0];
    leptons[3]=bosons[1]->children()[1];
  } else {
    leptons[2]=bosons[1]->children()[1];
    leptons[3]=bosons[1]->children()[0];
  }
  // For debugging; check that entries are what mcfm says they should be. 
//   if(leptons[0]->id()!=14||leptons[1]->id()!=-13||
//      leptons[2]->id()!=12||leptons[3]->id()!=-12) {
//     cout << "\n\n\n\n\n";
//     cout << "leptons[0] " << *leptons[0] << endl;
//     cout << "leptons[1] " << *leptons[1] << endl;
//     cout << "leptons[2] " << *leptons[2] << endl;
//     cout << "leptons[3] " << *leptons[3] << endl;
//   }

  AnalysisHandler::analyze( event, ieve, loop, state);
 
}

void VV_ME_Analysis::analyze(const tPVector & particles) {

  vector<Lorentz5Momentum> p;
  for(int i=0;i<=3;i++) p.push_back(leptons[i]->momentum());
  if(emitted) p.push_back(emitted->momentum());
  assert(p.size()==4||p.size()==5);

  int offset(-3);

  Lorentz5Momentum p34, p56, p3456;
  p34    = p[offset+3]+p[offset+4];
  p56    = p[offset+5]+p[offset+6];
  p3456  = p34+p56;
  Lorentz5Momentum p34rest(0.*GeV,0.*GeV,0.*GeV,p34.m(),p34.m());
  Lorentz5Momentum p3_Vrest, p4_Vrest;
  p3_Vrest = boostx(p[offset+3], p34, p34rest);
  p4_Vrest = boostx(p[offset+4], p34, p34rest);

  // cos(theta_p3) in the p34 rest frame:
  cos3_h.addWeighted(p3_Vrest.cosTheta(),1.);
  // cos(theta_p4) in the p34 rest frame:
  cos4_h.addWeighted(p4_Vrest.cosTheta(),1.);
  // Scalar sum of pts of p3 & p4:
  HT34_h.addWeighted((p[offset+3].perp()+p[offset+4].perp())/GeV,1.);
  // p3:
  eta3_h.addWeighted(p[offset+3].eta(),1.);
  pt3_h.addWeighted(p[offset+3].perp()/GeV,1.);
  // p4:
  eta4_h.addWeighted(p[offset+4].eta(),1.);
  pt4_h.addWeighted(p[offset+4].perp()/GeV,1.);
  // First vector boson:
  eta34_h.addWeighted(p34.eta(),1.);
  y34_h.addWeighted(p34.rapidity(),1.);
  pt34_h.addWeighted(p34.perp()/GeV,1.);
  m34wz_h.addWeighted(sqrt(p34.m2())/GeV,1.);
  // p5:
  eta5_h.addWeighted(p[offset+5].eta(),1.);
  pt5_h.addWeighted(p[offset+5].perp()/GeV,1.);
  // p6:
  eta6_h.addWeighted(p[offset+6].eta(),1.);
  pt6_h.addWeighted(p[offset+6].perp()/GeV,1.);
  // Second vector boson:
  m56_h.addWeighted(sqrt(p56.m2())/GeV,1.);
  eta56_h.addWeighted(p56.eta(),1.);
  y56_h.addWeighted(p56.rapidity(),1.);
  pt56_h.addWeighted(p56.perp()/GeV,1.);
  // Scalar sum of all lepton pts:
  HT3456_h.addWeighted((p[offset+3].perp()+p[offset+4].perp()
		      +p[offset+5].perp()+p[offset+6].perp()
		      )/GeV,1.);
  y3456_h.addWeighted(p3456.rapidity(),1.);
  m3456_h.addWeighted(sqrt(p3456.m2())/GeV,1.);
  // The theta Born variables:
  assert(abs(inbound.first->children()[0]->id())<7);
  assert(abs(inbound.second->children()[0]->id())<7);
  Lorentz5Momentum p3456rest(0.*GeV,0.*GeV,0.*GeV,p3456.m(),p3456.m());
  Lorentz5Momentum pplus(inbound.first->children()[0]->momentum()); 
  Lorentz5Momentum pminus(inbound.second->children()[0]->momentum());
  Lorentz5Momentum k(0.*GeV,0.*GeV,0.*GeV,0.*GeV,0.*GeV);
  pplus  = boostx(pplus , p3456, p3456rest);
  pminus = boostx(pminus, p3456, p3456rest);
  p34    = boostx(p34   , p3456, p3456rest);
  p56    = boostx(p56   , p3456, p3456rest);
  double cth1  (pplus.vect()*p34.vect()/pplus.vect().mag()/p34.vect().mag());
  cth1_h.addWeighted(cth1,1.);
  double sth1  (sqrt(1.-cth1*cth1));
  double cpsipr(0.);
  double spsipr(0.);
  double acth2 (0.);
  if(emitted) {
    k = emitted->momentum();
    k = boostx(k, p3456, p3456rest);
    cpsipr = pplus.vect()*k.vect()/pplus.vect().mag()/k.vect().mag(); 
    spsipr = sqrt(1.-cpsipr*cpsipr);
    acth2  = acos( ( k.vect()*p34.vect()/k.vect().mag()    /p34.vect().mag()
	           - cpsipr*cth1 
		   )/spsipr/sth1 );
  } else {
    acth2 = acos(p34.y()/p34.vect().mag()/sth1);
  }
  acth2_h.addWeighted(acth2,1.);

  AnalysisHandler::analyze(particles);
}

void VV_ME_Analysis::analyze(tPPtr part) {
  //find electron (and positron) with highest pt
  //check to see this is the same as the Z pt
  //   if( part->id() == ParticleID::eminus ) {
  //  cerr<< "id is: "<<part->id()<<" pt = "<< part->momentum().perp() / GeV<<" \n ";
  // }
}

void VV_ME_Analysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void VV_ME_Analysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<VV_ME_Analysis> VV_ME_Analysis::initVV_ME_Analysis;
// Definition of the static class description member.

void VV_ME_Analysis::Init() {

  static ClassDocumentation<VV_ME_Analysis> documentation
    ("There is no documentation for the VV_ME_Analysis class");

}

void VV_ME_Analysis::dofinish() {
  AnalysisHandler::dofinish();
  ofstream file;
  string fname = generator()->filename()+string("_")+name()+string(".top");

  // Normalise the histograms
  // cos(theta_p3) in the p34 rest frame:
  cos3_h.normaliseToCrossSection();
  cos3_h.prefactor(cos3_h.prefactor()*1.e6);
  // cos(theta_p4) in the p34 rest frame:
  cos4_h.normaliseToCrossSection();
  cos4_h.prefactor(cos4_h.prefactor()*1.e6);
  // Scalar sum of pts of p3 & p4:
  HT34_h.normaliseToCrossSection();
  HT34_h.prefactor(HT34_h.prefactor()*1.e6);
  // p3:
  eta3_h.normaliseToCrossSection();
  eta3_h.prefactor(eta3_h.prefactor()*1.e6);
  pt3_h.normaliseToCrossSection();
  pt3_h.prefactor(pt3_h.prefactor()*1.e6);
  // p4:
  eta4_h.normaliseToCrossSection();
  eta4_h.prefactor(eta4_h.prefactor()*1.e6);
  pt4_h.normaliseToCrossSection();
  pt4_h.prefactor(pt4_h.prefactor()*1.e6);
  // First vector boson:
  eta34_h.normaliseToCrossSection();
  eta34_h.prefactor(eta34_h.prefactor()*1.e6);
  y34_h.normaliseToCrossSection();
  y34_h.prefactor(y34_h.prefactor()*1.e6);
  pt34_h.normaliseToCrossSection();
  pt34_h.prefactor(pt34_h.prefactor()*1.e6);
  m34wz_h.normaliseToCrossSection();
  m34wz_h.prefactor(m34wz_h.prefactor()*1.e6);
  // p5:
  eta5_h.normaliseToCrossSection();
  eta5_h.prefactor(eta5_h.prefactor()*1.e6);
  pt5_h.normaliseToCrossSection();
  pt5_h.prefactor(pt5_h.prefactor()*1.e6);
  // p6:
  eta6_h.normaliseToCrossSection();
  eta6_h.prefactor(eta6_h.prefactor()*1.e6);
  pt6_h.normaliseToCrossSection();
  pt6_h.prefactor(pt6_h.prefactor()*1.e6);
  // Second vector boson:
  m56_h.normaliseToCrossSection();
  m56_h.prefactor(m56_h.prefactor()*1.e6);
  eta56_h.normaliseToCrossSection();
  eta56_h.prefactor(eta56_h.prefactor()*1.e6);
  y56_h.normaliseToCrossSection();
  y56_h.prefactor(y56_h.prefactor()*1.e6);
  pt56_h.normaliseToCrossSection();
  pt56_h.prefactor(pt56_h.prefactor()*1.e6);
  // Scalar sum of all lepton pts:
  HT3456_h.normaliseToCrossSection();
  HT3456_h.prefactor(HT3456_h.prefactor()*1.e6);
  y3456_h.normaliseToCrossSection();
  y3456_h.prefactor(y3456_h.prefactor()*1.e6);
  m3456_h.normaliseToCrossSection();
  m3456_h.prefactor(m3456_h.prefactor()*1.e6);
  // The theta Born variables:
  cth1_h.normaliseToCrossSection();
  cth1_h.prefactor(cth1_h.prefactor()*1.e6);
  acth2_h.normaliseToCrossSection();
  acth2_h.prefactor(acth2_h.prefactor()*1.e6);

  file.open(fname.c_str());
  using namespace HistogramOptions;

  // cos(theta_p3) in the p34 rest frame:
  // Lepton:
  cos3_h.topdrawOutput(file,Frame,"RED","cos3 distribution: all wgts");
  // cos(theta_p4) in the p34 rest frame:
  cos4_h.topdrawOutput(file,Frame,"RED","cos4 distribution: all wgts");
  // Scalar sum of pts of p3 & p4:
  HT34_h.topdrawOutput(file,Frame,"RED","HT34 distribution: all wgts");
  // p3:
  eta3_h.topdrawOutput(file,Frame,"RED","eta3 distribution: all wgts");
  pt3_h.topdrawOutput(file,Frame,"RED","pt3 distribution: all wgts");
  // p4:
  eta4_h.topdrawOutput(file,Frame,"RED","eta4 distribution: all wgts");
  pt4_h.topdrawOutput(file,Frame,"RED","pt4 distribution: all wgts");
  // First vector boson:
  eta34_h.topdrawOutput(file,Frame,"RED","eta34 distribution: all wgts");
  y34_h.topdrawOutput(file,Frame,"RED","y34 distribution: all wgts");
  pt34_h.topdrawOutput(file,Frame,"RED","pt34 low  mass: all wgts");
  m34wz_h.topdrawOutput(file,Frame,"RED","m34 low  mass: all wgts");
  // p5:
  eta5_h.topdrawOutput(file,Frame,"RED","eta5 distribution: all wgts");
  pt5_h.topdrawOutput(file,Frame,"RED","pt5 distribution: all wgts");
  // p6:
  eta6_h.topdrawOutput(file,Frame,"RED","eta6 distribution: all wgts");
  pt6_h.topdrawOutput(file,Frame,"RED","pt6 distribution: all wgts");
  // Second vector boson:
  m56_h.topdrawOutput(file,Frame,"RED","m56 distribution: all wgts");
  eta56_h.topdrawOutput(file,Frame,"RED","eta56 distribution: all wgts");
  y56_h.topdrawOutput(file,Frame,"RED","y56 distribution: all wgts");
  pt56_h.topdrawOutput(file,Frame,"RED","pt56 distribution: all wgts");
  // Scalar sum of all lepton pts:
  HT3456_h.topdrawOutput(file,Frame,"RED","HT3456 distribution: all wgts");
  y3456_h.topdrawOutput(file,Frame,"RED","y3456 distribution: all wgts");
  m3456_h.topdrawOutput(file,Frame,"RED","m3456 distribution: all wgts");
  // The theta Born variables:
  cth1_h.topdrawOutput(file,Frame,"RED","Cos(theta1) distribution: all wgts");
  acth2_h.topdrawOutput(file,Frame,"RED","Cos(theta2) distribution: all wgts");

  file.close();
}


Lorentz5Momentum VV_ME_Analysis::boostx(Lorentz5Momentum p_in, 
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
