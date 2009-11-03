// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VV_ME_Analysis class.
//

#include "VV_ME_Analysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

VV_ME_Analysis::VV_ME_Analysis() :
  _cos3_h(-1., 1.,40)    , _cos4_h(-1., 1.,40)  ,
  _cos5_h(-1., 1.,40)    , _cos6_h(-1., 1.,40)  ,
  _HT34_h(0.,500.,25)    ,
  _eta3_h(-6.,6.,100)    , _pt3_h(0.,200.,100)  ,
  _eta4_h(-6.,6.,100)    , _pt4_h(0.,200.,100)  ,
  _eta34_h(-6.,6.,100)   , _y34_h(-5.,5.,100)   ,
  _pt34_h(0.,500.,500)   , _pt34_mcfm_h(0.,500.,25),
  _m34wz_h(70.,100.,100) ,
  _eta5_h(-6.,6.,100)    , _pt5_h(0.,200.,100)  ,
  _eta6_h(-6.,6.,100)    , _pt6_h(0.,200.,100)  ,
  _eta56_h(-6.,6.,100)   , _y56_h(-5.,5.,100)   ,
  _pt56_h(0.,500.,500)   , _pt56_mcfm_h(0.,500.,25),
  _m56_h(75.0,100.0,100) ,
  _HT3456_h(0.,500.,25)  ,
  _y3456_h(-4.,4.,80)    , _m3456_h(0.,1000.,100),
  _th1_h(0.,Constants::pi,32), _th2_h(-Constants::pi,Constants::pi,32),
  _ptVV_h(0.,500.,500),_yVV_h(-4.,4.,80),
  _ptJet_h(0.,500.,500),
  _yJet_10_h(-10.0,10.0,200),_yJet_40_h(-10.0,10.0,200),
  _yJet_80_h(-10.0,10.0,200),
  _yJet_yV1_10_h(-10.0,10.0,200),_yJet_yV1_40_h(-10.0,10.0,200),
  _yJet_yV1_80_h(-10.0,10.0,200),
  _yJet_yV2_10_h(-10.0,10.0,200),_yJet_yV2_40_h(-10.0,10.0,200),
  _yJet_yV2_80_h(-10.0,10.0,200),
  _yJet_yVV_10_h(-10.0,10.0,200),_yJet_yVV_40_h(-10.0,10.0,200),
  _yJet_yVV_80_h(-10.0,10.0,200),
  _nJets_10_h(-0.5,100.5,101),_nJets_40_h(-0.5,100.5,101),
  _nJets_80_h(-0.5,100.5,101),
  _m34_mcatnlo_h(70.,100.,60),_m56_mcatnlo_h(70.,100.,60),
  _yJet2_10_h(-10.0,10.0,200),_yJet2_40_h(-10.0,10.0,200),
  _yJet2_80_h(-10.0,10.0,200),
  _yJet2_yVV_10_h(-10.0,10.0,200),_yJet2_yVV_40_h(-10.0,10.0,200),
  _yJet2_yVV_80_h(-10.0,10.0,200)
{}

void VV_ME_Analysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  // Store the initial state particles.
  inbound_  = event->incoming();

  // Store all of the final state particles as leptons.
  tPVector leptons;
  leptons = event->getFinalState();

  // Store all non-leptonic particles in partons.
  tPVector partons;
  for(unsigned int ix=0; ix<leptons.size();++ix)
    if(abs(leptons[ix]->id())>16||abs(leptons[ix]->id())<11)
      partons.push_back(leptons[ix]);

  // Delete all non-leptonic` particles from leptons.
  for(unsigned int ix=0; ix<leptons.size();++ix)
    if(abs(leptons[ix]->id())>16||abs(leptons[ix]->id())<11) { 
      leptons.erase(leptons.begin()+ix);
      ix--;
    }

  // Make sure that leptons has only 4 entries (2 bosons * 2-2-body decays).
  if(leptons.size()!=4) 
    throw Exception() << "VV_ME_Analysis::analyze\n"
		      << "Found " << leptons.size() << " leptons in final state!\n"
		      << Exception::eventerror;

  // Make sure all final-state leptons have parents.
  for (unsigned int i=0; i<leptons.size(); i++)
    assert(leptons[i]->parents().size()>0);

  // Go get the parents of the leptons W+ W- / W+ Z / W- Z / Z Z
  tPVector bosons;
  tPPtr aBoson;
  for(unsigned int ix=0;ix<leptons.size();ix++) {
    aBoson = leptons[ix]->parents()[0];
    while(leptons[ix]->id()==aBoson->id()&&leptons[ix]->parents().size()>0) 
      aBoson = aBoson->parents()[0];
    bosons.push_back(aBoson);
  }
  // Get rid of duplicate entries in bosons:
  tPVector::iterator pit;
  tPVector::iterator qit;
  for(pit=bosons.begin();pit!=bosons.end();++pit)
    for(qit=pit;qit!=bosons.end();++qit)
      if((*qit)->id()==(*pit)->id()) {
	bosons.erase(qit);
	break;
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
//     cout << "\n\n\n\n\n";
//     cout << "leptons[0] " << *leptons[0] << endl;
//     cout << "leptons[1] " << *leptons[1] << endl;
//     cout << "leptons[2] " << *leptons[2] << endl;
//     cout << "leptons[3] " << *leptons[3] << endl;

  // Get the jets ordered by their Pt (largest Pt is first).
  vector<fastjet::PseudoJet> particlesToCluster;
  for(unsigned int jx=0; jx<partons.size();jx++) {
    fastjet::PseudoJet p(partons[jx]->momentum().x()/GeV,
			 partons[jx]->momentum().y()/GeV,
			 partons[jx]->momentum().z()/GeV,
			 partons[jx]->momentum().e()/GeV);
    p.set_user_index(jx);
    particlesToCluster.push_back(p);
  }
  fastjet::RecombinationScheme recombinationScheme = fastjet::E_scheme;
  fastjet::Strategy            strategy            = fastjet::Best;
  double R(0.7);
  fastjet::JetDefinition       jetDefinition(fastjet::kt_algorithm,
					     R,
					     recombinationScheme,
					     strategy);
  fastjet::ClusterSequence fastjetEvent(particlesToCluster,jetDefinition);
  vector<fastjet::PseudoJet> inclusiveJets = fastjetEvent.inclusive_jets();
  inclusiveJets = fastjet::sorted_by_pt(inclusiveJets);

  epsilon_ = 1.e-10;

  vector<Lorentz5Momentum> p;
  for(int i=0;i<=3;i++) p.push_back(leptons[i]->momentum());

  int offset(-3);

  Lorentz5Momentum p34, p56, p3456;
  p34    = p[offset+3]+p[offset+4];
  p56    = p[offset+5]+p[offset+6];
  p3456  = p34+p56;
  Lorentz5Momentum p34rest(0.*GeV,0.*GeV,0.*GeV,p34.m(),p34.m());
  Lorentz5Momentum p56rest(0.*GeV,0.*GeV,0.*GeV,p56.m(),p56.m());
  Lorentz5Momentum p3_Vrest, p4_Vrest, p5_Vrest, p6_Vrest;
  p3_Vrest = boostx(p[offset+3], p34, p34rest);
  p4_Vrest = boostx(p[offset+4], p34, p34rest);
  p5_Vrest = boostx(p[offset+5], p56, p56rest);
  p6_Vrest = boostx(p[offset+6], p56, p56rest);

  // cos(theta_p3) in the p34 rest frame:
  _cos3_h.addWeighted(p3_Vrest.cosTheta(),1.);
  // cos(theta_p4) in the p34 rest frame:
  _cos4_h.addWeighted(p4_Vrest.cosTheta(),1.);
  // cos(theta_p5) in the p56 rest frame:
  _cos5_h.addWeighted(p5_Vrest.cosTheta(),1.);
  // cos(theta_p6) in the p56 rest frame:
  _cos6_h.addWeighted(p6_Vrest.cosTheta(),1.);
  // Scalar sum of pts of p3 & p4:
  _HT34_h.addWeighted((p[offset+3].perp()+p[offset+4].perp())/GeV,1.);
  // p3:
  _eta3_h.addWeighted(p[offset+3].eta(),1.);
  _pt3_h.addWeighted(p[offset+3].perp()/GeV,1.);
  // p4:
  _eta4_h.addWeighted(p[offset+4].eta(),1.);
  _pt4_h.addWeighted(p[offset+4].perp()/GeV,1.);
  // First vector boson:
  _eta34_h.addWeighted(p34.eta(),1.);
  _y34_h.addWeighted(p34.rapidity(),1.);
  _pt34_h.addWeighted(p34.perp()/GeV,1.);
  _pt34_mcfm_h.addWeighted(p34.perp()/GeV,1.);
  _m34wz_h.addWeighted(sqrt(p34.m2())/GeV,1.);
  // p5:
  _eta5_h.addWeighted(p[offset+5].eta(),1.);
  _pt5_h.addWeighted(p[offset+5].perp()/GeV,1.);
  // p6:
  _eta6_h.addWeighted(p[offset+6].eta(),1.);
  _pt6_h.addWeighted(p[offset+6].perp()/GeV,1.);
  // Second vector boson:
  _eta56_h.addWeighted(p56.eta(),1.);
  _y56_h.addWeighted(p56.rapidity(),1.);
  _pt56_h.addWeighted(p56.perp()/GeV,1.);
  _pt56_mcfm_h.addWeighted(p56.perp()/GeV,1.);
  _m56_h.addWeighted(sqrt(p56.m2())/GeV,1.);
  // The two vector bosons together:
  _ptVV_h.addWeighted((p34+p56).perp()/GeV,1.);
  _yVV_h.addWeighted((p34+p56).rapidity(),1.);
  // The hardest jet stuff - pT followed by rapidity:
  if(inclusiveJets.size()>0) {
    fastjet::PseudoJet hardestJet(inclusiveJets[0]);
    _ptJet_h.addWeighted(hardestJet.perp(),1.);
    if(hardestJet.perp()>10.) {
      _yJet_10_h.addWeighted(hardestJet.rapidity(),1.);
      _yJet_yV1_10_h.addWeighted(hardestJet.rapidity()-p34.rapidity(),1.);
      _yJet_yV2_10_h.addWeighted(hardestJet.rapidity()-p56.rapidity(),1.);
      _yJet_yVV_10_h.addWeighted(hardestJet.rapidity()-(p34+p56).rapidity(),1.);
    }
    if(hardestJet.perp()>40.) {
      _yJet_40_h.addWeighted(hardestJet.rapidity(),1.);
      _yJet_yV1_40_h.addWeighted(hardestJet.rapidity()-p34.rapidity(),1.);
      _yJet_yV2_40_h.addWeighted(hardestJet.rapidity()-p56.rapidity(),1.);
      _yJet_yVV_40_h.addWeighted(hardestJet.rapidity()-(p34+p56).rapidity(),1.);
    } 
    if(hardestJet.perp()>80.) {
      _yJet_80_h.addWeighted(hardestJet.rapidity(),1.);
      _yJet_yV1_80_h.addWeighted(hardestJet.rapidity()-p34.rapidity(),1.);
      _yJet_yV2_80_h.addWeighted(hardestJet.rapidity()-p56.rapidity(),1.);
      _yJet_yVV_80_h.addWeighted(hardestJet.rapidity()-(p34+p56).rapidity(),1.);
    }
    if(inclusiveJets.size()>1) {
      fastjet::PseudoJet secondhardestJet(inclusiveJets[1]);
      if(secondhardestJet.perp()>10.) {
	_yJet2_10_h.addWeighted(secondhardestJet.rapidity(),1.);
	_yJet2_yVV_10_h.addWeighted(secondhardestJet.rapidity()-(p34+p56).rapidity(),1.);
      }
      if(secondhardestJet.perp()>40.) {
	_yJet2_40_h.addWeighted(secondhardestJet.rapidity(),1.);
	_yJet2_yVV_40_h.addWeighted(secondhardestJet.rapidity()-(p34+p56).rapidity(),1.);
      } 
      if(secondhardestJet.perp()>80.) {
	_yJet2_80_h.addWeighted(secondhardestJet.rapidity(),1.);
	_yJet2_yVV_80_h.addWeighted(secondhardestJet.rapidity()-(p34+p56).rapidity(),1.);
      }
    }
  }
  unsigned int n10(0);
  unsigned int n40(0);
  unsigned int n80(0);
  for(unsigned int ix=0; ix<inclusiveJets.size(); ++ix) {
    if(inclusiveJets[ix].perp()>=10.) n10++;
    if(inclusiveJets[ix].perp()>=40.) n40++;
    if(inclusiveJets[ix].perp()>=80.) n80++;
  }
  _nJets_10_h.addWeighted(n10+0.001,1.);
  _nJets_40_h.addWeighted(n40+0.001,1.);
  _nJets_80_h.addWeighted(n80+0.001,1.);
  // Masses of the two vector bosons:
  _m34_mcatnlo_h.addWeighted(p34.m()/GeV,1.);
  _m56_mcatnlo_h.addWeighted(p56.m()/GeV,1.);
  // Scalar sum of all lepton pts:
  _HT3456_h.addWeighted((p[offset+3].perp()+p[offset+4].perp()
		        +p[offset+5].perp()+p[offset+6].perp()
		       )/GeV,1.);
  _y3456_h.addWeighted(p3456.rapidity(),1.);
  _m3456_h.addWeighted(sqrt(p3456.m2())/GeV,1.);
  // The theta Born variables:
  assert(abs(inbound_.first->children()[0]->id())<7);
  assert(abs(inbound_.second->children()[0]->id())<7);
  Lorentz5Momentum pplus(inbound_.first->children()[0]->momentum()); 
  Lorentz5Momentum pminus(inbound_.second->children()[0]->momentum());
  Lorentz5Momentum p3456rest(0.*GeV,0.*GeV,0.*GeV,p3456.m(),p3456.m());
  // Boost everything to the diboson rest frame:
  pplus  = boostx(pplus , p3456, p3456rest);
  pminus = boostx(pminus, p3456, p3456rest);
  p34    = boostx(p34   , p3456, p3456rest);
  p56    = boostx(p56   , p3456, p3456rest);
  double sinthpplus (pplus.perp() /pplus.vect().mag() );
  double sinthpminus(pminus.perp()/pminus.vect().mag());
  if(partons.size()==0&&(sinthpplus>epsilon_||sinthpminus>epsilon_)) 
    throw Exception() << "VV_ME_Analysis::analyze\n"
		      << "Found 0 QCD particles in final state but\n" 
		      << "sin(theta_pplus)  = " << sinthpplus  << "\n"
		      << "sin(theta_pminus) = " << sinthpminus << "\n"
		      << Exception::warning;
  // Get the rotation that puts pplus on the z-axis (zrot):
  if(sinthpplus>epsilon_&&partons.size()>0) {
    LorentzRotation zrot(hwurot(pplus,1.,0.));
    // Now rotate everything using this matrix
    pplus  *= zrot;
    pminus *= zrot;
    p34    *= zrot;
    p56    *= zrot;
  }

  if(sinthpminus>epsilon_&&partons.size()>0) {
    // Get the rotation that puts the transverse bit 
    // of pminus on the +y axis (xyrot):
    LorentzRotation xyrot(hwurot(pplus, pminus.x()/pminus.perp(),
				        pminus.y()/pminus.perp()));
    
    // Now rotate everything using this matrix
    pplus  *= xyrot;
    pminus *= xyrot;
    p34    *= xyrot;
    p56    *= xyrot;
  }
  // Now get the Born variables
  _th1_h.addWeighted(acos(p34.z()/p34.vect().mag()),1.);
  _th2_h.addWeighted(atan2(p34.x(),p34.y()),1.);

}

void VV_ME_Analysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.

}

void VV_ME_Analysis::analyze(tPPtr) {}

// ClassDescription<VV_ME_Analysis> VV_ME_Analysis::initVV_ME_Analysis;
// // Definition of the static class description member.

NoPIOClassDescription<VV_ME_Analysis> VV_ME_Analysis::initVV_ME_Analysis;
// Definition of the static class description member.

void VV_ME_Analysis::Init() {

  static ClassDocumentation<VV_ME_Analysis> documentation
    ("There is no documentation for the VV_ME_Analysis class");

}

void VV_ME_Analysis::dofinish() {
  AnalysisHandler::dofinish();

  using namespace HistogramOptions;

  ofstream mcfm_file;
  string mcfm_filename 
    = generator()->filename()+string("_")+name()+string("_mcfm.top");
  mcfm_file.open(mcfm_filename.c_str());

  ofstream mcatnlo_file;
  string mcatnlo_filename 
    = generator()->filename()+string("_")+name()+string("_mcatnlo.top");
  mcatnlo_file.open(mcatnlo_filename.c_str());

  // cos(theta_p3) in the p34 rest frame:
  _cos3_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos3 distribution Unit Area");
  _cos3_h.normaliseToCrossSection();
  _cos3_h.prefactor(_cos3_h.prefactor()*1.e6);
  _cos3_h.topdrawOutput(mcfm_file,Frame,"RED","cos3 distribution: all wgts");
  _cos3_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos3 distribution cross");
  // cos(theta_p4) in the p34 rest frame:
  _cos4_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos4 distribution Unit Area");
  _cos4_h.normaliseToCrossSection();
  _cos4_h.prefactor(_cos4_h.prefactor()*1.e6);
  _cos4_h.topdrawOutput(mcfm_file,Frame,"RED","cos4 distribution: all wgts");
  _cos4_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos4 distribution cross");
  // cos(theta_p5) in the p56 rest frame:
  _cos5_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos5 distribution Unit Area");
  _cos5_h.normaliseToCrossSection();
  _cos5_h.prefactor(_cos5_h.prefactor()*1.e6);
  _cos5_h.topdrawOutput(mcfm_file,Frame,"RED","cos5 distribution: all wgts");
  _cos5_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos5 distribution cross");
  // cos(theta_p6) in the p56 rest frame:
  _cos6_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos6 distribution Unit Area");
  _cos6_h.normaliseToCrossSection();
  _cos6_h.prefactor(_cos6_h.prefactor()*1.e6);
  _cos6_h.topdrawOutput(mcfm_file,Frame,"RED","cos6 distribution: all wgts");
  _cos6_h.topdrawOutput(mcatnlo_file,Frame,"RED","cos6 distribution cross");
  // Scalar sum of pts of p3 & p4:
  _HT34_h.normaliseToCrossSection();
  _HT34_h.prefactor(_HT34_h.prefactor()*1.e6);
  _HT34_h.topdrawOutput(mcfm_file,Frame,"RED","HT34 distribution: all wgts");

  // p3:
  _eta3_h.normaliseToCrossSection();
  _eta3_h.prefactor(_eta3_h.prefactor()*1.e6);
  _eta3_h.topdrawOutput(mcfm_file,Frame,"RED","eta3 distribution: all wgts");
  _pt3_h.normaliseToCrossSection();
  _pt3_h.prefactor(_pt3_h.prefactor()*1.e6);
  _pt3_h.topdrawOutput(mcfm_file,Frame,"RED","pt3 distribution: all wgts");
  // p4:
  _eta4_h.normaliseToCrossSection();
  _eta4_h.prefactor(_eta4_h.prefactor()*1.e6);
  _eta4_h.topdrawOutput(mcfm_file,Frame,"RED","eta4 distribution: all wgts");
  _pt4_h.normaliseToCrossSection();
  _pt4_h.prefactor(_pt4_h.prefactor()*1.e6);
  _pt4_h.topdrawOutput(mcfm_file,Frame,"RED","pt4 distribution: all wgts");

  // First vector boson:
  _eta34_h.topdrawOutput(mcatnlo_file,Frame,"RED","eta34 Unit Area");
  _eta34_h.normaliseToCrossSection();
  _eta34_h.prefactor(_eta34_h.prefactor()*1.e6);
  _eta34_h.topdrawOutput(mcatnlo_file,Frame,"RED","eta34 cross");

  _eta34_h.normaliseToCrossSection();
  _eta34_h.prefactor(_eta34_h.prefactor()*1.e6);
  _eta34_h.topdrawOutput(mcfm_file,Frame,"RED","eta34 distribution: all wgts");

  _y34_h.topdrawOutput(mcatnlo_file,Frame,"RED","y34 Unit Area");
  _y34_h.normaliseToCrossSection();
  _y34_h.prefactor(_y34_h.prefactor()*1.e6);
  _y34_h.topdrawOutput(mcatnlo_file,Frame,"RED","y34 cross");

  _y34_h.normaliseToCrossSection();
  _y34_h.prefactor(_y34_h.prefactor()*1.e6);
  _y34_h.topdrawOutput(mcfm_file,Frame,"RED","y34 distribution: all wgts");

  _pt34_h.topdrawOutput(mcatnlo_file,Frame,"RED","V1 Pt Unit Area");
  _pt34_h.normaliseToCrossSection();
  _pt34_h.prefactor(_pt34_h.prefactor()*1.e6);
  _pt34_h.topdrawOutput(mcatnlo_file,Frame,"RED","V1 Pt cross");

  _pt34_mcfm_h.normaliseToCrossSection();
  _pt34_mcfm_h.prefactor(_pt34_mcfm_h.prefactor()*1.e6);
  _pt34_mcfm_h.topdrawOutput(mcfm_file,Frame,"RED","pt34 low  mass: all wgts");

  _m34wz_h.normaliseToCrossSection();
  _m34wz_h.prefactor(_m34wz_h.prefactor()*1.e6);
  _m34wz_h.topdrawOutput(mcfm_file,Frame,"RED","m34 low  mass: all wgts");

  // p5:
  _eta5_h.normaliseToCrossSection();
  _eta5_h.prefactor(_eta5_h.prefactor()*1.e6);
  _eta5_h.topdrawOutput(mcfm_file,Frame,"RED","eta5 distribution: all wgts");
  _pt5_h.normaliseToCrossSection();
  _pt5_h.prefactor(_pt5_h.prefactor()*1.e6);
  _pt5_h.topdrawOutput(mcfm_file,Frame,"RED","pt5 distribution: all wgts");
  // p6:
  _eta6_h.normaliseToCrossSection();
  _eta6_h.prefactor(_eta6_h.prefactor()*1.e6);
  _eta6_h.topdrawOutput(mcfm_file,Frame,"RED","eta6 distribution: all wgts");
  _pt6_h.normaliseToCrossSection();
  _pt6_h.prefactor(_pt6_h.prefactor()*1.e6);
  _pt6_h.topdrawOutput(mcfm_file,Frame,"RED","pt6 distribution: all wgts");

  // Second vector boson:
  _eta56_h.topdrawOutput(mcatnlo_file,Frame,"RED","eta56 Unit Area");
  _eta56_h.normaliseToCrossSection();
  _eta56_h.prefactor(_eta56_h.prefactor()*1.e6);
  _eta56_h.topdrawOutput(mcatnlo_file,Frame,"RED","eta56 cross");

  _y56_h.topdrawOutput(mcatnlo_file,Frame,"RED","y56 Unit Area");
  _y56_h.normaliseToCrossSection();
  _y56_h.prefactor(_y56_h.prefactor()*1.e6);
  _y56_h.topdrawOutput(mcatnlo_file,Frame,"RED","y56 cross");

  _pt56_h.topdrawOutput(mcatnlo_file,Frame,"RED","V2 Pt Unit Area");
  _pt56_h.normaliseToCrossSection();
  _pt56_h.prefactor(_pt56_h.prefactor()*1.e6);
  _pt56_h.topdrawOutput(mcatnlo_file,Frame,"RED","V2 Pt cross");

  _m56_h.normaliseToCrossSection();
  _m56_h.prefactor(_m56_h.prefactor()*1.e6);
  _m56_h.topdrawOutput(mcfm_file,Frame,"RED","m56 distribution: all wgts");

  _eta56_h.normaliseToCrossSection();
  _eta56_h.prefactor(_eta56_h.prefactor()*1.e6);
  _eta56_h.topdrawOutput(mcfm_file,Frame,"RED","eta56 distribution: all wgts");

  _y56_h.normaliseToCrossSection();
  _y56_h.prefactor(_y56_h.prefactor()*1.e6);
  _y56_h.topdrawOutput(mcfm_file,Frame,"RED","y56 distribution: all wgts");

  _pt56_mcfm_h.normaliseToCrossSection();
  _pt56_mcfm_h.prefactor(_pt56_mcfm_h.prefactor()*1.e6);
  _pt56_mcfm_h.topdrawOutput(mcfm_file,Frame,"RED","pt56 low  mass: all wgts");

  // The two vector bosons together:
  _ptVV_h.topdrawOutput(mcatnlo_file,Frame,"RED","VV Pt Unit Area");
  _ptVV_h.normaliseToCrossSection();
  _ptVV_h.prefactor(_ptVV_h.prefactor()*1.e6);
  _ptVV_h.topdrawOutput(mcatnlo_file,Frame,"RED","VV Pt cross");
  _yVV_h.topdrawOutput(mcatnlo_file,Frame,"RED","VV rapidity Unit Area");
  _yVV_h.normaliseToCrossSection();
  _yVV_h.prefactor(_yVV_h.prefactor()*1.e6);
  _yVV_h.topdrawOutput(mcatnlo_file,Frame,"RED","VV rapidity cross");

  // The hardest jet stuff - pT followed by rapidity:
  _ptJet_h.topdrawOutput(mcatnlo_file,Frame,"RED","Hardest Jet pT Unit Area");
  _ptJet_h.normaliseToCrossSection();
  _ptJet_h.prefactor(_ptJet_h.prefactor()*1.e6);
  _ptJet_h.topdrawOutput(mcatnlo_file,Frame,"RED","Hardest Jet pT cross");
  _yJet_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>10 Unit Area");
  _yJet_10_h.normaliseToCrossSection();
  _yJet_10_h.prefactor(_yJet_10_h.prefactor()*1.e6);
  _yJet_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>10 cross");
  _yJet_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>40 Unit Area");
  _yJet_40_h.normaliseToCrossSection();
  _yJet_40_h.prefactor(_yJet_40_h.prefactor()*1.e6);
  _yJet_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>40 cross");
  _yJet_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>80 Unit Area");
  _yJet_80_h.normaliseToCrossSection();
  _yJet_80_h.prefactor(_yJet_80_h.prefactor()*1.e6);
  _yJet_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet pT>80 cross");

  // The hardest jet rapidity gap w.r.t V1:
  _yJet_yV1_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>10 Unit Area");
  _yJet_yV1_10_h.normaliseToCrossSection();
  _yJet_yV1_10_h.prefactor(_yJet_yV1_10_h.prefactor()*1.e6);
  _yJet_yV1_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>10 cross");
  _yJet_yV1_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>40 Unit Area");
  _yJet_yV1_40_h.normaliseToCrossSection();
  _yJet_yV1_40_h.prefactor(_yJet_yV1_40_h.prefactor()*1.e6);
  _yJet_yV1_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>40 cross");
  _yJet_yV1_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>80 Unit Area");
  _yJet_yV1_80_h.normaliseToCrossSection();
  _yJet_yV1_80_h.prefactor(_yJet_yV1_80_h.prefactor()*1.e6);
  _yJet_yV1_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV1 pT>80 cross");

  // The hardest jet rapidity gap w.r.t V2:
  _yJet_yV2_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>10 Unit Area");
  _yJet_yV2_10_h.normaliseToCrossSection();
  _yJet_yV2_10_h.prefactor(_yJet_yV2_10_h.prefactor()*1.e6);
  _yJet_yV2_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>10 cross");
  _yJet_yV2_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>40 Unit Area");
  _yJet_yV2_40_h.normaliseToCrossSection();
  _yJet_yV2_40_h.prefactor(_yJet_yV2_40_h.prefactor()*1.e6);
  _yJet_yV2_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>40 cross");
  _yJet_yV2_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>80 Unit Area");
  _yJet_yV2_80_h.normaliseToCrossSection();
  _yJet_yV2_80_h.prefactor(_yJet_yV2_80_h.prefactor()*1.e6);
  _yJet_yV2_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yV2 pT>80 cross");

  // The hardest jet rapidity gap w.r.t VV:
  _yJet_yVV_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>10 Unit Area");
  _yJet_yVV_10_h.normaliseToCrossSection();
  _yJet_yVV_10_h.prefactor(_yJet_yVV_10_h.prefactor()*1.e6);
  _yJet_yVV_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>10 cross");
  _yJet_yVV_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>40 Unit Area");
  _yJet_yVV_40_h.normaliseToCrossSection();
  _yJet_yVV_40_h.prefactor(_yJet_yVV_40_h.prefactor()*1.e6);
  _yJet_yVV_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>40 cross");
  _yJet_yVV_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>80 Unit Area");
  _yJet_yVV_80_h.normaliseToCrossSection();
  _yJet_yVV_80_h.prefactor(_yJet_yVV_80_h.prefactor()*1.e6);
  _yJet_yVV_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet - yVV pT>80 cross");

  // Jet multiplicities:
  _nJets_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>10 Unit Area");
  _nJets_10_h.normaliseToCrossSection();
  _nJets_10_h.prefactor(_nJets_10_h.prefactor()*1.e6);
  _nJets_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>10 cross");
  _nJets_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>40 Unit Area");
  _nJets_40_h.normaliseToCrossSection();
  _nJets_40_h.prefactor(_nJets_40_h.prefactor()*1.e6);
  _nJets_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>40 cross");
  _nJets_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>80 Unit Area");
  _nJets_80_h.normaliseToCrossSection();
  _nJets_80_h.prefactor(_nJets_80_h.prefactor()*1.e6);
  _nJets_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","N Jets pT>80 cross");

  // Masses of the two vector bosons:
  _m34_mcatnlo_h.topdrawOutput(mcatnlo_file,Frame,"RED","m34 Unit Area");
  _m34_mcatnlo_h.normaliseToCrossSection();
  _m34_mcatnlo_h.prefactor(_m34_mcatnlo_h.prefactor()*1.e6);
  _m34_mcatnlo_h.topdrawOutput(mcatnlo_file,Frame,"RED","m34 cross");
  _m56_mcatnlo_h.topdrawOutput(mcatnlo_file,Frame,"RED","m56 Unit Area");
  _m56_mcatnlo_h.normaliseToCrossSection();
  _m56_mcatnlo_h.prefactor(_m56_mcatnlo_h.prefactor()*1.e6);
  _m56_mcatnlo_h.topdrawOutput(mcatnlo_file,Frame,"RED","m56 cross");

  // The second hardest jet rapidity:
  _yJet2_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>10 Unit Area");
  _yJet2_10_h.normaliseToCrossSection();
  _yJet2_10_h.prefactor(_yJet2_10_h.prefactor()*1.e6);
  _yJet2_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>10 cross");
  _yJet2_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>40 Unit Area");
  _yJet2_40_h.normaliseToCrossSection();
  _yJet2_40_h.prefactor(_yJet2_40_h.prefactor()*1.e6);
  _yJet2_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>40 cross");
  _yJet2_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>80 Unit Area");
  _yJet2_80_h.normaliseToCrossSection();
  _yJet2_80_h.prefactor(_yJet2_80_h.prefactor()*1.e6);
  _yJet2_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","y Jet2 pT>80 cross");

  // The second hardest jet rapidity gap w.r.t VV:
  _yJet2_yVV_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>10 Unit Area");
  _yJet2_yVV_10_h.normaliseToCrossSection();
  _yJet2_yVV_10_h.prefactor(_yJet2_yVV_10_h.prefactor()*1.e6);
  _yJet2_yVV_10_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>10 cross");
  _yJet2_yVV_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>40 Unit Area");
  _yJet2_yVV_40_h.normaliseToCrossSection();
  _yJet2_yVV_40_h.prefactor(_yJet2_yVV_40_h.prefactor()*1.e6);
  _yJet2_yVV_40_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>40 cross");
  _yJet2_yVV_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>80 Unit Area");
  _yJet2_yVV_80_h.normaliseToCrossSection();
  _yJet2_yVV_80_h.prefactor(_yJet2_yVV_80_h.prefactor()*1.e6);
  _yJet2_yVV_80_h.topdrawOutput(mcatnlo_file,Frame,"RED","yJet2 - yVV pT>80 cross");

  // Scalar sum of all lepton pts:
  _HT3456_h.normaliseToCrossSection();
  _HT3456_h.prefactor(_HT3456_h.prefactor()*1.e6);
  _HT3456_h.topdrawOutput(mcfm_file,Frame,"RED","HT3456 distribution: all wgts");
  _y3456_h.normaliseToCrossSection();
  _y3456_h.prefactor(_y3456_h.prefactor()*1.e6);
  _y3456_h.topdrawOutput(mcfm_file,Frame,"RED","y3456 distribution: all wgts");
  _m3456_h.normaliseToCrossSection();
  _m3456_h.prefactor(_m3456_h.prefactor()*1.e6);
  _m3456_h.topdrawOutput(mcfm_file,Frame,"RED","m3456 distribution: all wgts");

  // The theta Born variables:
  _th1_h.normaliseToCrossSection();
  _th1_h.prefactor(_th1_h.prefactor()*1.e6);
  _th1_h.topdrawOutput(mcfm_file,Frame,"RED","theta1 distribution: all wgts");
  _th2_h.normaliseToCrossSection();
  _th2_h.prefactor(_th2_h.prefactor()*1.e6);
  _th2_h.topdrawOutput(mcfm_file,Frame,"RED","theta2 distribution: all wgts");

  mcfm_file.close();
  mcatnlo_file.close();
}


Lorentz5Momentum VV_ME_Analysis::boostx(Lorentz5Momentum p_in, 
		  		        Lorentz5Momentum pt, 
				        Lorentz5Momentum ptt){
  // Boost input vector p_in to output vector p_out using the same
  // transformation as required to boost massive vector pt to ptt
  Lorentz5Momentum p_tmp,p_out;
  ThreeVector<double> beta;
  Energy  mass, bdotp;
  double  gam;

  if (pt.m2()<0.*GeV2)
    throw Exception() << "pt.m2()<0. in boostx, pt.m2() = "
		      << pt.m2()/GeV2 
		      << Exception::runerror;

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

LorentzRotation VV_ME_Analysis::hwurot(Lorentz5Momentum p,double cp,double sp) {
  //-----------------------------------------------------------------------
  // r is a  rotation matrix to get from vector p the z-axis, followed by
  // a rotation by phi about the z-axis, where cp = cos(phi), sp = sin(phi)
  //-----------------------------------------------------------------------

  LorentzRotation rxy;
  LorentzRotation ryz;

  LorentzRotation raz;
  p.z()/GeV > 0 ? raz.setRotateZ( atan2(cp,sp)) 
                : raz.setRotateZ(-atan2(cp,sp));

  if(p.perp()/p.vect().mag()>=epsilon_) {
    rxy.setRotateZ(atan2(p.x(),p.y()));
    p *= rxy;
    ryz.setRotateX(acos(p.z()/p.vect().mag()));
    p *= ryz;
    return raz*ryz*rxy;
  } 
  else {
    return raz;
  }
  
}
