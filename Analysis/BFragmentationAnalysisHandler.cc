// -*- C++ -*-
//
// BFragmentationAnalysisHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BFragmentationAnalysisHandler class.
//

#include "BFragmentationAnalysisHandler.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig/Utilities/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void BFragmentationAnalysisHandler::analyze(tEventPtr event, long,
					    int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  double weight = event->weight();
  ///////////////////////////
  // Hadron Level Analysis //
  ///////////////////////////
  // extract the weakly decaying B hadrons using set to avoid double counting
  set<PPtr> allParticles;
  event->select(inserter(allParticles),WeakBHadronSelector());
  // convert to vector
  tPVector particles(allParticles.begin(),allParticles.end());
  // numerator
  _emax = 0.5*generator()->maximumCMEnergy();   
  analyze(particles,weight); 

  ////////////////////////////////////////
  // Parton Level Analysis (e+e-->bbar) //
  ////////////////////////////////////////
  // Get all the particles from the perturbative bit of the events. Find the
  // b and bbar coming straight out of the Z/photon (b_orig, bbar_orig and 
  // ZGamma respectively). Then go and find the b and bbar that have their 
  // maximal off-shellness i.e. just before the first gluon is emitted. 
  // Finally go off down the b quark lines to find the b's just before they
  // hadronize (b_end and bbar_end respectively).
  ParticleSet pert=event->primaryCollision()->step(1)->all();
  analyze_bquarks(pert,weight);
}

void BFragmentationAnalysisHandler::analyze(tPPtr part, double weight) {
  _fragBxE ->addWeighted(part->momentum().e()/_emax, weight);
  _fragBxEa->addWeighted(part->momentum().e()/_emax, weight);
}

void BFragmentationAnalysisHandler::analyze_bquarks(ParticleSet pert, 
						    double weight) {
  ParticleSet::const_iterator pit;
  PPtr b_orig,bbar_orig,ZGamma,b_start,bbar_start,b_end,bbar_end;
  // First go through all the particles looking for b's coming out of Z/gamma's:
  for(pit=pert.begin();pit!=pert.end();++pit) {
    PPtr bline = *pit;
    if(abs((*pit)->id())==5) {
      while(abs(bline->parents()[0]->id())==5) bline=bline->parents()[0];  
    }
    if(bline->id()== 5&&(bline->parents()[0]->id()!=21)) b_orig    = bline;
    if(bline->id()==-5&&(bline->parents()[0]->id()!=21)) bbar_orig = bline;
  }
  if(!b_orig) return;
  if(!bbar_orig) return;
  // Note down the Z/Photon that decays to the b & bbar:
  if(b_orig->parents()[0]==bbar_orig->parents()[0]) 
      ZGamma = b_orig->parents()[0];
  PPtr root_b[] = {b_orig,bbar_orig};
  // Now go and look for the b & bbar just before the first gluon is emitted:
  for(int ix=0;ix<=1;ix++) {
      while(root_b[ix]->momentum().m()<=
	    getParticleData(ParticleID::b)->mass()+1.e-8*GeV) {
	  for(unsigned int jx=0;jx<root_b[ix]->children().size();jx++) 
	      if(root_b[ix]->id()==root_b[ix]->children()[jx]->id()) 
		  root_b[ix]=root_b[ix]->children()[jx]; 
      }
  }
  b_start    = root_b[0];
  bbar_start = root_b[1];
  // Now go and find the b and bbar quarks at the end of the shower before 
  // they turn into hadrons.
  for(unsigned int ix=0;ix<=1;ix++) {
      while(root_b[ix]->momentum().m()>=
	    getParticleData(ParticleID::b)->constituentMass()+1.e-8*GeV) {
	  for(unsigned int jx=0;jx<root_b[ix]->children().size();jx++) 
	      if(root_b[ix]->id()==root_b[ix]->children()[jx]->id()) 
		  root_b[ix]=root_b[ix]->children()[jx]; 
      }
  }
  b_end    = root_b[0];
  bbar_end = root_b[1];
  // Fill the energy fraction histograms with that of the b quarks.
  _fragbquarkxE      ->
    addWeighted( b_end->momentum().e()/_emax   ,weight);
  _fragbquarkxE      ->
    addWeighted( bbar_end->momentum().e()/_emax,weight);
  _fragbquarkjetmass ->
    addWeighted( b_start->momentum().m()/GeV   ,weight);
  _fragbquarkjetmass ->
    addWeighted( bbar_start->momentum().m()/GeV,weight);
}

NoPIOClassDescription<BFragmentationAnalysisHandler> 
BFragmentationAnalysisHandler::initBFragmentationAnalysisHandler;
// Definition of the static class description member.

void BFragmentationAnalysisHandler::Init() {

  static ClassDocumentation<BFragmentationAnalysisHandler> documentation
    ("The BFragmentationAnalysisHandler class performs analysis"
     " of the B fragmentation function",
     "The B fragmentation function analysis uses data from \\cite{Heister:2001jg,Abe:2002iq}.",
     "  %\\cite{Heister:2001jg}\n"
     "\\bibitem{Heister:2001jg}\n"
     "  A.~Heister {\\it et al.}  [ALEPH Collaboration],\n"
     "  %``Study of the fragmentation of b quarks into B mesons at the Z peak,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 512}, 30 (2001)\n"
     "  [arXiv:hep-ex/0106051].\n"
     "  %%CITATION = PHLTA,B512,30;%%\n"
     "%\\cite{Abe:2002iq}\n"
     "\\bibitem{Abe:2002iq}\n"
     "  K.~Abe {\\it et al.}  [SLD Collaboration],\n"
     "  %``Measurement of the b-quark fragmentation function in Z0 decays,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 65}, 092006 (2002)\n"
     "  [Erratum-ibid.\\  D {\\bf 66}, 079905 (2002)]\n"
     "  [arXiv:hep-ex/0202031].\n"
     "  %%CITATION = PHRVA,D65,092006;%%\n"
     );

}

void BFragmentationAnalysisHandler::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  // output the histograms
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  double chisq,minfrac=0.05;
  unsigned int npoint;
  _fragBxE->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD b hadron fragmentation "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _fragBxE->topdrawOutput(output,Frame|Errorbars,
			  "RED",
			  "B Hadron fragmentation function compared to SLD data",
			  "                                                    ",
			  "1/SdS/dx0B1",
			  "  G G   X X",
			  "x0B1",
			  " X X");
  _fragBxEa->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for ALEPH b hadron fragmentation "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _fragBxEa->topdrawOutput(output,Frame|Errorbars,
			   "RED",
			   "B Hadron framgentation function compared to ALEPH data",
			   "                                                      ",
			   "1/SdS/dx0B1",
			   "  G G   X X",
			   "x0B1",
			   " X X");
  _fragbquarkxE->topdrawOutput(output,Frame|Smooth,
		       "RED",
		       "b quark fragmentation function",
		       "                              ",
		       "1/SdS/dx0B1",
		       "  G G   X X",
		       "x0B1",
		       " X X");
  _fragbquarkjetmass->topdrawOutput(output,Frame|Smooth,
		       "RED",
		       "b quark jet mass",
		       "                                                    ",
		       "1/SdS/dm0J1223",
		       "  G G   X XX X",
		       "m0J1223",
		       " X XX X");
  output.close();
}

void BFragmentationAnalysisHandler::doinitrun() {
  AnalysisHandler::doinitrun();
  // SLD binning
  double BxEbins[] = {0.00, 0.04, 0.08, 0.12, 0.16, 
		      0.20, 0.24, 0.28, 0.32, 0.36, 
		      0.40, 0.44, 0.48, 0.52, 0.56, 
		      0.60, 0.64, 0.68, 0.72, 0.76, 
		      0.80, 0.84, 0.88, 0.92, 0.96, 
		      1.0};
  double BxEdata[] = {0.000,0.000,0.000,0.116,0.198,
		      0.247,0.264,0.308,0.370,0.426,
		      0.501,0.577,0.685,0.833,1.055,
		      1.311,1.667,2.080,2.566,2.934,
		      3.104,2.856,1.954,0.841,0.108};
  double BxEerror[]= {0.000,0.000,0.000,0.030,0.037,
		      0.030,0.029,0.032,0.033,0.034,
		      0.039,0.041,0.042,0.053,0.074,
		      0.089,0.088,0.084,0.116,0.178,
		      0.235,0.179,0.162,0.215,0.062};
  vector<double> bins(BxEbins,BxEbins+26), data(BxEdata,BxEdata+25),
    error(BxEerror,BxEerror+25);
  _fragBxE  = new_ptr(Histogram(bins,data,error));
  // ALEPH binning
  double BxEabins[] = {0.0   ,0.1   ,0.25  ,0.35  ,0.45  , 
		       0.55  ,0.6   ,0.65  ,0.7   ,0.725 , 
		       0.75  ,0.775 ,0.8   ,0.825 ,0.85  , 
		       0.875 ,0.9   ,0.925 ,0.95  ,0.975 , 
		       1.};
  double BxEadata[] = {0.0000,0.1193,0.2810,0.4510,0.7410,
		       1.0180,1.2760,1.7020,2.1080,2.3520,
		       2.5360,2.7960,2.9840,3.1000,2.9080,
		       2.6440,2.0880,1.3480,0.4840,0.0400};
  double BxEaerrora[] = {0.0000,0.0487,0.0470,0.0390,0.0590,
			 0.0660,0.0640,0.0660,0.0760,0.0840,
			 0.0920,0.1040,0.1080,0.1000,0.0880,
			 0.0880,0.1160,0.1240,0.0760,0.0120};
  double BxEaerrorb[] = {0.0000,0.0573,0.0350,0.0430,0.0660,
			 0.0680,0.0640,0.0740,0.0960,0.1120,
			 0.1240,0.1360,0.1320,0.1240,0.1040,
			 0.1360,0.1880,0.1880,0.1000,0.0200};
  double BxEaerror[20];
  for(unsigned int ix=0;ix<20;++ix){BxEaerror[ix]=sqrt(sqr(BxEaerrora[ix])+
						       sqr(BxEaerrorb[ix]));}
  bins  = vector<double>(BxEabins,BxEabins+21);
  data  = vector<double>(BxEadata,BxEadata+20);
  error = vector<double>(BxEaerror,BxEaerror+20);
  _fragBxEa = new_ptr(Histogram(bins,data,error));
  _fragbquarkxE = new_ptr(Histogram(0.,1.0,100));
  _fragbquarkjetmass = new_ptr(Histogram(0.,90.0,90));
}
