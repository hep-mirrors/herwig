 // -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PP2HAnalysis class.
//

#include "PP2HAnalysis.h"
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

Histogram pt_30(  0., 30.,50);
Histogram pt_200( 0.,200.,50);
Histogram pt_400( 0.,400.,50);
Histogram Yjet_10(   -5.,5.,50);
Histogram YjetYH_10( -5.,5.,50);
Histogram Yjet_40(   -5.,5.,50);
Histogram YjetYH_40( -5.,5.,50);
Histogram Yjet_80(   -5.,5.,50);
Histogram YjetYH_80( -5.,5.,50);
Histogram Yjet_10A(  -5.,5.,50);
Histogram YjetYH_10A(-5.,5.,50);
Histogram Yjet_40A(  -5.,5.,50);
Histogram YjetYH_40A(-5.,5.,50);
Histogram Yjet_80A(  -5.,5.,50);
Histogram YjetYH_80A(-5.,5.,50);
Histogram Njets_10(0.,10.,10);
Histogram Njets_40(0.,10.,10);
Histogram Njets_80(0.,10.,10);
Histogram log_y23(-14.,-4.,100);
Histogram log_y34(-14.,-4.,100);
Histogram log_y45(-14.,-4.,100);
Histogram cosjet(-1.,1.,50);
Histogram cosjet_lt_10(-1.,1.,50);
Histogram cosjet_lt_30(-1.,1.,50);
Histogram cosjet_lt_100(-1.,1.,50);
Histogram cosjet_gt_10(-1.,1.,50);
Histogram cosjet_gt_30(-1.,1.,50);
Histogram cosjet_gt_100(-1.,1.,50);
Histogram cosjet_w(-0.5,0.5,50);
Histogram cosjet_lt_10_w(-0.5,0.5,50);
Histogram cosjet_lt_30_w(-0.5,0.5,50);
Histogram cosjet_lt_100_w(-0.5,0.5,50);
Histogram cosjet_gt_10_w(-0.5,0.5,50);
Histogram cosjet_gt_30_w(-0.5,0.5,50);
Histogram cosjet_gt_100_w(-0.5,0.5,50);
Histogram cos1st(-1.,1.,50);
Histogram cos1st_lt_10(-1.,1.,50);
Histogram cos1st_lt_30(-1.,1.,50);
Histogram cos1st_lt_100(-1.,1.,50);
Histogram cos1st_gt_10(-1.,1.,50);
Histogram cos1st_gt_30(-1.,1.,50);
Histogram cos1st_gt_100(-1.,1.,50);
Histogram cos1st_w(-0.5,0.5,50);
Histogram cos1st_lt_10_w(-0.5,0.5,50);
Histogram cos1st_lt_30_w(-0.5,0.5,50);
Histogram cos1st_lt_100_w(-0.5,0.5,50);
Histogram cos1st_gt_10_w(-0.5,0.5,50);
Histogram cos1st_gt_30_w(-0.5,0.5,50);
Histogram cos1st_gt_100_w(-0.5,0.5,50);

void PP2HAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void PP2HAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<PP2HAnalysis> PP2HAnalysis::initPP2HAnalysis;
// Definition of the static class description member.

void PP2HAnalysis::Init() {

  static ClassDocumentation<PP2HAnalysis> documentation
    ("There is no documentation for the PP2HAnalysis class");

}

void PP2HAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {

  // Store all of the final state particles as particles.
  particles=event->getFinalState();
  // Extract the emitted parton
  tPVector::iterator pit;
  int Higgses(0);
  for(pit=particles.begin();pit!=particles.end();++pit) {
    if((*pit)->id()==25) { 
       Higgs = *pit;
       particles.erase(pit); // Delete the Higgs so it is not clustered
       Higgses++;
    }
  }
  if(Higgses>1||Higgses==0)
	throw Exception() << "PP2HAnalysis::analyze"  
			  << "\nMore / less than one Higgs in the final state"
			  << "\nNumberof Higgses is " << Higgses 
			  << Exception::warning;

  //-----------------------//
  //                       //
  // pT Spectrum analysis  //
  //                       //
  //-----------------------//
  Energy pTh = Higgs->momentum().perp();
  pt_30.addWeighted( pTh/GeV,1.);
  pt_200.addWeighted(pTh/GeV,1.);
  pt_400.addWeighted(pTh/GeV,1.);

  //-----------------------//
  //                       //
  // Jet Rapidity analysis //
  //                       //
  //-----------------------//
  _kint.clearMap();
  KtJet::KtEvent ev = KtJet::KtEvent(_kint.convert(particles),4,2,1,0.7);
  // Get the two jets ordered by their Pt (largest Pt is first).
  vector<KtJet::KtLorentzVector> ktjets = ev.getJetsPt();
  if( ktjets.size()>0
      && Higgs
      && abs(ktjets[0].t())!=abs(ktjets[0].z())
      && abs(Higgs->momentum().t())!=abs(Higgs->momentum().z()) ) {
    double Yjet    = ktjets[0].rapidity();
    double Yjet_YH = ktjets[0].rapidity()-Higgs->momentum().rapidity();
    if(ktjets[0].perp()>=10000.) {
      Yjet_10.addWeighted(Yjet,1.);
      YjetYH_10.addWeighted(Yjet_YH,1.);
    }
    if(ktjets[0].perp()>=40000.) {
      Yjet_40.addWeighted(Yjet,1.);
      YjetYH_40.addWeighted(Yjet_YH,1.);
    }
    if(ktjets[0].perp()>=80000.) {
      Yjet_80.addWeighted(Yjet,1.);
      YjetYH_80.addWeighted(Yjet_YH,1.);
    }
    unsigned int last(ktjets.size()-1);
    if(ktjets[last]>=10000.&&Higgs->momentum().perp()>=10000.*MeV) {
      Yjet_10A.addWeighted(Yjet,1.);
      YjetYH_10A.addWeighted(Yjet_YH,1.);
    }
    if(ktjets[last]>=40000.&&Higgs->momentum().perp()>=40000.*MeV) {
      Yjet_40A.addWeighted(Yjet,1.);
      YjetYH_40A.addWeighted(Yjet_YH,1.);
    }
    if(ktjets[last]>=80000.&&Higgs->momentum().perp()>=80000.*MeV) {
      Yjet_80A.addWeighted(Yjet,1.);
      YjetYH_80A.addWeighted(Yjet_YH,1.);
    }
  } else {
    if(!Higgs) 
      generator()->log() << "PP2HAnalysis::analyze\n"
			 << "Didn't find a Higgs in final state!\n";
    if(ktjets.size()==0)
      generator()->log()   << "PP2HAnalysis::analyze\n"
			 << "Didn't find any jets!\n";
    if(abs(ktjets[0].t())!=abs(ktjets[0].z()))
      generator()->log()   << "PP2HAnalysis::analyze\n"
			 << "Infinite rapidity leading jet.\n";
  }

  //----------------------//
  //                      //
  //    Soft radiation    //
  //       analysis       //
  //                      //
  //----------------------//
  if( ktjets.size()>0 && Higgs ) {

    // Find the gluons which fuse to give the Higgs:
    tPPtr parent_P;
    tPPtr parent_M;
    if(Higgs->parents().size()==1 && 
       Higgs->parents()[0]->parents().size()==2 && 
       Higgs->parents()[0]->parents()[0]->id()==21 &&
       Higgs->parents()[0]->parents()[1]->id()==21) {
      parent_P = Higgs->parents()[0]->parents()[0];
      parent_M = Higgs->parents()[0]->parents()[1];
    } else if(Higgs->parents().size()==2 &&
	      Higgs->parents()[0]->id()==21 && 
	      Higgs->parents()[1]->id()==21) {
      parent_P = Higgs->parents()[0];
      parent_M = Higgs->parents()[1];
    } else {
      throw Exception() << "PP2HAnalysis::analyze"  
			<< "\nCan't find Higgs parent gluons"
			<< Exception::warning;
    }
    if(parent_P->momentum().z()<0.*GeV) swap(parent_P,parent_M);

    // Find the boost to the COM of the colliding gluons.
    Lorentz5Momentum p_P_gg = parent_P->momentum();
    Lorentz5Momentum p_M_gg = parent_M->momentum();
    Lorentz5Momentum pjet_gg = _kint.convert(ktjets[0]);
    Boost boost_gg = (p_P_gg+p_M_gg).findBoostToCM();
    // Boost the two gluons to their COM with them aligned along Z. 
    p_P_gg = p_P_gg.boost(boost_gg);
    p_M_gg = p_M_gg.boost(boost_gg);
    // Also boost the leading jet to this frame.
    pjet_gg = pjet_gg.boost(boost_gg);

    // Find the first emission (with the widest angle).
    tPPtr emitter_P = parent_P;
    tPPtr ISchild_P;
    tPPtr emitted_P;
    tPPtr emitter_M = parent_M;
    tPPtr ISchild_M;
    tPPtr emitted_M;
    do {
      ISchild_P = emitter_P;
      emitter_P = emitter_P->parents()[0];
    }
    while(emitter_P->children().size()==1);
    emitted_P=emitter_P->children()[0]==ISchild_P ? emitter_P->children()[1] :
                                                    emitter_P->children()[0] ;
    do {
      ISchild_M = emitter_M;
      emitter_M = emitter_M->parents()[0];
    }
    while(emitter_M->children().size()==1);
    emitted_M=emitter_M->children()[0]==ISchild_M ? emitter_M->children()[1] :
                                                    emitter_M->children()[0] ;

    if(abs(emitter_P->id())==2212||abs(emitter_M->id())==2212) 
      throw Exception() << "PP2HAnalysis::analyze"  
			<< "\nSelected a proton as a parent of emission."
			<< Exception::warning;
      
    Lorentz5Momentum p1st_P_gg = emitted_P->momentum();
    Lorentz5Momentum p1st_M_gg = emitted_M->momentum();
    p1st_P_gg = p1st_P_gg.boost(boost_gg);
    p1st_M_gg = p1st_M_gg.boost(boost_gg);

    double c1st_P = p1st_P_gg.cosTheta();
    double c1st_M = p1st_M_gg.cosTheta();
    double c1st   = abs(c1st_P) <= abs(c1st_M) ? c1st_P : c1st_M;
    if(abs(emitter_P->id())==2212) c1st = c1st_M;
    if(abs(emitter_M->id())==2212) c1st = c1st_P;
    double cjet   = pjet_gg.cosTheta();
    cos1st.addWeighted(c1st,1.);
    cosjet.addWeighted(cjet,1.);
    cos1st_w.addWeighted(c1st,1.);
    cosjet_w.addWeighted(cjet,1.);
    if(ktjets[0].perp()<=10000.) {
      cos1st_lt_10.addWeighted(c1st,1.);
      cosjet_lt_10.addWeighted(cjet,1.);
      cos1st_lt_10_w.addWeighted(c1st,1.);
      cosjet_lt_10_w.addWeighted(cjet,1.);
    }
    if(ktjets[0].perp()<=30000.) {
      cos1st_lt_30.addWeighted(c1st,1.);
      cosjet_lt_30.addWeighted(cjet,1.);
      cos1st_lt_30_w.addWeighted(c1st,1.);
      cosjet_lt_30_w.addWeighted(cjet,1.);
    }
    if(ktjets[0].perp()<=100000.) {
      cos1st_lt_100.addWeighted(c1st,1.);
      cosjet_lt_100.addWeighted(cjet,1.);
      cos1st_lt_100_w.addWeighted(c1st,1.);
      cosjet_lt_100_w.addWeighted(cjet,1.);
    }
    if(ktjets[0].perp()>=10000.) {
      cos1st_gt_10.addWeighted(c1st,1.);
      cosjet_gt_10.addWeighted(cjet,1.);
      cos1st_gt_10_w.addWeighted(c1st,1.);
      cosjet_gt_10_w.addWeighted(cjet,1.);
    }
    if(ktjets[0].perp()>=30000.) {
      cos1st_gt_30.addWeighted(c1st,1.);
      cosjet_gt_30.addWeighted(cjet,1.);
      cos1st_gt_30_w.addWeighted(c1st,1.);
      cosjet_gt_30_w.addWeighted(cjet,1.);
    }
    if(ktjets[0].perp()>=100000.) {
      cos1st_gt_100.addWeighted(c1st,1.);
      cosjet_gt_100.addWeighted(cjet,1.);
      cos1st_gt_100_w.addWeighted(c1st,1.);
      cosjet_gt_100_w.addWeighted(cjet,1.);
    }
  }

  //----------------------//
  //                      //
  //   Jet multiplicity   //
  //       analysis       //
  //                      //
  //----------------------//
  unsigned int n10(0);
  unsigned int n40(0);
  unsigned int n80(0);
  for(unsigned int ix=0; ix<ktjets.size(); ++ix) {
    if(ktjets[ix].perp()>=10000.) n10++;
    if(ktjets[ix].perp()>=40000.) n40++;
    if(ktjets[ix].perp()>=80000.) n80++;
  }
  Njets_10.addWeighted(n10+0.001,1);
  Njets_40.addWeighted(n40+0.001,1);
  Njets_80.addWeighted(n80+0.001,1);

  //----------------------//
  //                      //
  //   Y merge analysis   //
  //                      //
  //----------------------//
  if(ktjets.size()>2)
    log_y23.addWeighted(log(ev.getYMerge(2))/log(10.),1); 
  if(ktjets.size()>3)
    log_y34.addWeighted(log(ev.getYMerge(3))/log(10.),1); 
  if(ktjets.size()>4)
    log_y45.addWeighted(log(ev.getYMerge(4))/log(10.),1); 

}

void PP2HAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename()+string("-")+name()+string(".top");
  ofstream outfile(fname.c_str());

  using namespace HistogramOptions;
  pt_30.normaliseToCrossSection();
  pt_30.prefactor(pt_30.prefactor()*1.e3);
  pt_30.topdrawOutput(outfile,Frame,"RED",
		      "Higgs boson p0T1",
		      "             X X",
		      "dS/dp0T1 (pb)",
		      " G   X X     ");
  pt_200.normaliseToCrossSection();
  pt_200.prefactor(pt_200.prefactor()*1.e3);
  pt_200.topdrawOutput(outfile,Frame|Ylog,"RED",
		       "Higgs boson p0T1",
		       "             X X",
		       "dS/dp0T1 (pb)",
		       " G   X X     ");
  pt_400.normaliseToCrossSection();
  pt_400.prefactor(pt_400.prefactor()*1.e3);
  pt_400.topdrawOutput(outfile,Frame|Ylog,"RED",
		       "Higgs boson p0T1",
		       "             X X",
		       "dS/dp0T1 (pb)",
		       " G   X X     ");
  Yjet_10.normaliseToCrossSection();
  Yjet_10.prefactor(Yjet_10.prefactor()*1.e3);
  Yjet_10.topdrawOutput(outfile,Frame|Ylog,"RED",
			"y0jet1 (p0T,jet1>10 GeV)",
			" X   X   X     X        ",
			"dS/dy0jet1 (pb)",
			" G   X   X     ");
  Yjet_40.normaliseToCrossSection();
  Yjet_40.prefactor(Yjet_40.prefactor()*1.e3);
  Yjet_40.topdrawOutput(outfile,Frame|Ylog,"RED",
			"y0jet1 (p0T,jet1>40 GeV)",
			" X   X   X     X        ",
			"dS/dy0jet1 (pb)",
			" G   X   X     ");
  Yjet_80.normaliseToCrossSection();
  Yjet_80.prefactor(Yjet_80.prefactor()*1.e3);
  Yjet_80.topdrawOutput(outfile,Frame|Ylog,"RED",
			"y0jet1 (p0T,jet1>80 GeV)",
			" X   X   X     X        ",
			"dS/dy0jet1 (pb)",
			" G   X   X     ");
  YjetYH_10.normaliseToCrossSection();
  YjetYH_10.prefactor(YjetYH_10.prefactor()*1.e3);
  YjetYH_10.topdrawOutput(outfile,Frame|Ylog,"RED",
			  "y0jet1-y0H1 (p0T,jet1>10 GeV)",
			  " X   X  X X   X     X        ",
			  "dS/d(y0jet1-y0H1) (pb)",
			  " G    X   X  X X      ");
  YjetYH_40.normaliseToCrossSection();
  YjetYH_40.prefactor(YjetYH_40.prefactor()*1.e3);
  YjetYH_40.topdrawOutput(outfile,Frame|Ylog,"RED",
			  "y0jet1-y0H1 (p0T,jet1>40 GeV)",
			  " X   X  X X   X     X        ",
			  "dS/d(y0jet1-y0H1) (pb)",
			  " G    X   X  X X      ");
  YjetYH_80.normaliseToCrossSection();
  YjetYH_80.prefactor(YjetYH_80.prefactor()*1.e3);
  YjetYH_80.topdrawOutput(outfile,Frame|Ylog,"RED",
			  "y0jet1-y0H1 (p0T,jet1>80 GeV)",
			  " X   X  X X   X     X        ",
			  "dS/d(y0jet1-y0H1) (pb)",
			  " G    X   X  X X      ");
  Yjet_10A.normaliseToCrossSection();
  Yjet_10A.prefactor(Yjet_10A.prefactor()*1.e3);
  Yjet_10A.topdrawOutput(outfile,Frame|Ylog,"RED",
			 "y0jet1 (p0T,all1>10 GeV)",
			 " X   X   X     X        ",
			 "dS/dy0jet1 (pb)",
			 " G   X   X     ");
  Yjet_40A.normaliseToCrossSection();
  Yjet_40A.prefactor(Yjet_40A.prefactor()*1.e3);
  Yjet_40A.topdrawOutput(outfile,Frame|Ylog,"RED",
			 "y0jet1 (p0T,all1>40 GeV)",
			 " X   X   X     X        ",
			 "dS/dy0jet1 (pb)",
			 " G   X   X     ");
  Yjet_80A.normaliseToCrossSection();
  Yjet_80A.prefactor(Yjet_80A.prefactor()*1.e3);
  Yjet_80A.topdrawOutput(outfile,Frame|Ylog,"RED",
			 "y0jet1 (p0T,all1>80 GeV)",
			 " X   X   X     X        ",
			 "dS/dy0jet1 (pb)",
			 " G   X   X     ");
  YjetYH_10A.normaliseToCrossSection();
  YjetYH_10A.prefactor(YjetYH_10A.prefactor()*1.e3);
  YjetYH_10A.topdrawOutput(outfile,Frame|Ylog,"RED",
			   "y0jet1-y0H1 (p0T,all1>10 GeV)",
			   " X   X  X X   X     X        ",
			   "dS/d(y0jet1-y0H1) (pb)",
			   " G    X   X  X X      ");
  YjetYH_40A.normaliseToCrossSection();
  YjetYH_40A.prefactor(YjetYH_40A.prefactor()*1.e3);
  YjetYH_40A.topdrawOutput(outfile,Frame|Ylog,"RED",
			   "y0jet1-y0H1 (p0T,all1>40 GeV)",
			   " X   X  X X   X     X        ",
			   "dS/d(y0jet1-y0H1) (pb)",
			   " G    X   X  X X      ");
  YjetYH_80A.normaliseToCrossSection();
  YjetYH_80A.prefactor(YjetYH_80A.prefactor()*1.e3);
  YjetYH_80A.topdrawOutput(outfile,Frame|Ylog,"RED",
			   "y0jet1-y0H1 (p0T,all1>80 GeV)",
			   " X   X  X X   X     X        ",
			   "dS/d(y0jet1-y0H1) (pb)",
			   " G    X   X  X X      ");
  Njets_10.normaliseToCrossSection();
  Njets_10.prefactor(Njets_10.prefactor()*1.e3);
  Njets_10.topdrawOutput(outfile,Frame,"RED",
			 "Jet multiplicity (p0T,jets1>10 GeV)",
			 "                   X      X        ",
			 "dS/dn0jets1 (pb)",
			 " G   X    X     ");
  Njets_40.normaliseToCrossSection();
  Njets_40.prefactor(Njets_40.prefactor()*1.e3);
  Njets_40.topdrawOutput(outfile,Frame,"RED",
			 "Jet multiplicity (p0T,jets1>40 GeV)",
			 "                   X      X        ",
			 "dS/dn0jets1 (pb)",
			 " G   X    X     ");
  Njets_80.normaliseToCrossSection();
  Njets_80.prefactor(Njets_80.prefactor()*1.e3);
  Njets_80.topdrawOutput(outfile,Frame,"RED",
			 "Jet multiplicity (p0T,jets1>80 GeV)",
			 "                   X      X        ",
			 "dS/dn0jets1 (pb)",
			 " G   X    X     ");

  log_y23.normaliseToCrossSection();
  log_y23.prefactor(log_y23.prefactor()*1.e3);
  log_y23.topdrawOutput(outfile,Frame|Ylog,"RED",
			"Log0101(y0231)",
			"   X  X  X  X ",
			"dS/dLog0101(y0231)",
			" G     X  X  X  X ");
  log_y34.normaliseToCrossSection();
  log_y34.prefactor(log_y34.prefactor()*1.e3);
  log_y34.topdrawOutput(outfile,Frame|Ylog,"RED",
			"Log0101(y0341)",
			"   X  X  X  X ",
			"dS/dLog0101(y0341)",
			" G     X  X  X  X ");
  log_y45.normaliseToCrossSection();
  log_y45.prefactor(log_y45.prefactor()*1.e3);
  log_y45.topdrawOutput(outfile,Frame|Ylog,"RED",
			"Log0101(y0451)",
			"   X  X  X  X ",
			"dS/dLog0101(y0451)",
			" G     X  X  X  X ");
  cosjet.normaliseToCrossSection();
  cosjet.prefactor(cosjet.prefactor()*1.e3);
  cosjet.topdrawOutput(outfile,Frame,"RED",
		       "CosQ0jet1",
		       "   GX   X",
		       "dS/dcosQ0jet1 (pb)",
		       " G     GX   X     ");
  cosjet_lt_10.normaliseToCrossSection();
  cosjet_lt_10.prefactor(cosjet_lt_10.prefactor()*1.e3);
  cosjet_lt_10.topdrawOutput(outfile,Frame,"RED",
			     "CosQ0jet1 (p0t,jet1<10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ0jet1 (pb)",
			     " G     GX   X     ");
  cosjet_lt_30.normaliseToCrossSection();
  cosjet_lt_30.prefactor(cosjet_lt_30.prefactor()*1.e3);
  cosjet_lt_30.topdrawOutput(outfile,Frame,"RED",
			     "CosQ0jet1 (p0t,jet1<30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ0jet1 (pb)",
			     " G     GX   X     ");
  cosjet_lt_100.normaliseToCrossSection();
  cosjet_lt_100.prefactor(cosjet_lt_100.prefactor()*1.e3);
  cosjet_lt_100.topdrawOutput(outfile,Frame,"RED",
			      "CosQ0jet1 (p0t,jet1<100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ0jet1 (pb)",
			      " G     GX   X     ");
  cosjet_gt_10.normaliseToCrossSection();
  cosjet_gt_10.prefactor(cosjet_gt_10.prefactor()*1.e3);
  cosjet_gt_10.topdrawOutput(outfile,Frame,"RED",
			     "CosQ0jet1 (p0t,jet1>10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ0jet1 (pb)",
			     " G     GX   X     ");
  cosjet_gt_30.normaliseToCrossSection();
  cosjet_gt_30.prefactor(cosjet_gt_30.prefactor()*1.e3);
  cosjet_gt_30.topdrawOutput(outfile,Frame,"RED",
			     "CosQ0jet1 (p0t,jet1>30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ0jet1 (pb)",
			     " G     GX   X     ");
  cosjet_gt_100.normaliseToCrossSection();
  cosjet_gt_100.prefactor(cosjet_gt_100.prefactor()*1.e3);
  cosjet_gt_100.topdrawOutput(outfile,Frame,"RED",
			      "CosQ0jet1 (p0t,jet1>100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ0jet1 (pb)",
			      " G     GX   X     ");
  cosjet_lt_10_w.normaliseToCrossSection();
  cosjet_lt_10_w.prefactor(cosjet_lt_10_w.prefactor()*1.e3);
  cosjet_lt_10_w.topdrawOutput(outfile,Frame,"RED",
			       "CosQ0jet1 (p0t,jet1<10 GeV)",
			       "   GX   X   X     X        ",
			       "dS/dcosQ0jet1 (pb)",
			       " G     GX   X     ");
  cosjet_lt_30_w.normaliseToCrossSection();
  cosjet_lt_30_w.prefactor(cosjet_lt_30_w.prefactor()*1.e3);
  cosjet_lt_30_w.topdrawOutput(outfile,Frame,"RED",
			       "CosQ0jet1 (p0t,jet1<30 GeV)",
			       "   GX   X   X     X        ",
			       "dS/dcosQ0jet1 (pb)",
			       " G     GX   X     ");
  cosjet_lt_100_w.normaliseToCrossSection();
  cosjet_lt_100_w.prefactor(cosjet_lt_100_w.prefactor()*1.e3);
  cosjet_lt_100_w.topdrawOutput(outfile,Frame,"RED",
				"CosQ0jet1 (p0t,jet1<100 GeV)",
				"   GX   X   X     X        ",
				"dS/dcosQ0jet1 (pb)",
				" G     GX   X     ");
  cosjet_gt_10_w.normaliseToCrossSection();
  cosjet_gt_10_w.prefactor(cosjet_gt_10_w.prefactor()*1.e3);
  cosjet_gt_10_w.topdrawOutput(outfile,Frame,"RED",
			       "CosQ0jet1 (p0t,jet1>10 GeV)",
			       "   GX   X   X     X        ",
			       "dS/dcosQ0jet1 (pb)",
			       " G     GX   X     ");
  cosjet_gt_30_w.normaliseToCrossSection();
  cosjet_gt_30_w.prefactor(cosjet_gt_30_w.prefactor()*1.e3);
  cosjet_gt_30_w.topdrawOutput(outfile,Frame,"RED",
			       "CosQ0jet1 (p0t,jet1>30 GeV)",
			       "   GX   X   X     X        ",
			       "dS/dcosQ0jet1 (pb)",
			       " G     GX   X     ");
  cosjet_gt_100_w.normaliseToCrossSection();
  cosjet_gt_100_w.prefactor(cosjet_gt_100_w.prefactor()*1.e3);
  cosjet_gt_100_w.topdrawOutput(outfile,Frame,"RED",
				"CosQ0jet1 (p0t,jet1>100 GeV)",
				"   GX   X   X     X        ",
				"dS/dcosQ0jet1 (pb)",
				" G     GX   X     ");
  cos1st.normaliseToCrossSection();
  cos1st.prefactor(cos1st.prefactor()*1.e3);
  cos1st.topdrawOutput(outfile,Frame,"RED",
		       "CosQ01st1",
		       "   GX   X",
		       "dS/dcosQ01st1 (pb)",
		       " G     GX   X     ");
  cos1st_lt_10.normaliseToCrossSection();
  cos1st_lt_10.prefactor(cos1st_lt_10.prefactor()*1.e3);
  cos1st_lt_10.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1<10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_lt_30.normaliseToCrossSection();
  cos1st_lt_30.prefactor(cos1st_lt_30.prefactor()*1.e3);
  cos1st_lt_30.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1<30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_lt_100.normaliseToCrossSection();
  cos1st_lt_100.prefactor(cos1st_lt_100.prefactor()*1.e3);
  cos1st_lt_100.topdrawOutput(outfile,Frame,"RED",
			      "CosQ01st1 (p0T,jet1<100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ01st1 (pb)",
			      " G     GX   X     ");
  cos1st_gt_10.normaliseToCrossSection();
  cos1st_gt_10.prefactor(cos1st_gt_10.prefactor()*1.e3);
  cos1st_gt_10.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1>10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_gt_30.normaliseToCrossSection();
  cos1st_gt_30.prefactor(cos1st_gt_30.prefactor()*1.e3);
  cos1st_gt_30.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1>30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_gt_100.normaliseToCrossSection();
  cos1st_gt_100.prefactor(cos1st_gt_100.prefactor()*1.e3);
  cos1st_gt_100.topdrawOutput(outfile,Frame,"RED",
			      "CosQ01st1 (p0T,jet1>100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ01st1 (pb)",
			      " G     GX   X     ");
  cos1st_lt_10_w.normaliseToCrossSection();
  cos1st_lt_10_w.prefactor(cos1st_lt_10_w.prefactor()*1.e3);
  cos1st_lt_10_w.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1<10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_lt_30_w.normaliseToCrossSection();
  cos1st_lt_30_w.prefactor(cos1st_lt_30_w.prefactor()*1.e3);
  cos1st_lt_30_w.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1<30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_lt_100_w.normaliseToCrossSection();
  cos1st_lt_100_w.prefactor(cos1st_lt_100_w.prefactor()*1.e3);
  cos1st_lt_100_w.topdrawOutput(outfile,Frame,"RED",
			      "CosQ01st1 (p0T,jet1<100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ01st1 (pb)",
			      " G     GX   X     ");
  cos1st_gt_10_w.normaliseToCrossSection();
  cos1st_gt_10_w.prefactor(cos1st_gt_10_w.prefactor()*1.e3);
  cos1st_gt_10_w.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1>10 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_gt_30_w.normaliseToCrossSection();
  cos1st_gt_30_w.prefactor(cos1st_gt_30_w.prefactor()*1.e3);
  cos1st_gt_30_w.topdrawOutput(outfile,Frame,"RED",
			     "CosQ01st1 (p0T,jet1>30 GeV)",
			     "   GX   X   X     X        ",
			     "dS/dcosQ01st1 (pb)",
			     " G     GX   X     ");
  cos1st_gt_100_w.normaliseToCrossSection();
  cos1st_gt_100_w.prefactor(cos1st_gt_100_w.prefactor()*1.e3);
  cos1st_gt_100_w.topdrawOutput(outfile,Frame,"RED",
			      "CosQ01st1 (p0T,jet1>100 GeV)",
			      "   GX   X   X     X        ",
			      "dS/dcosQ01st1 (pb)",
			      " G     GX   X     ");

  outfile.close();
}

