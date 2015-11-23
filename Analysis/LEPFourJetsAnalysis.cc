// -*- C++ -*-
//
// LEPFourJetsAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LepFourJetsAnalysis class.
//

#include "LEPFourJetsAnalysis.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

namespace {

  struct ChargedFinalState {
    static bool AllCollisions() { return false; }
    static bool AllSteps() { return false; }
    // ===
    // pick the last instance from the shower
    static bool FinalState() { return true; }
    static bool Intermediate() { return false; }
    // ===
    static bool Check(const Particle & p) { 
      return ParticleTraits<Particle>::iCharge(p);
    }
  };
  
}

void LEPFourJetsAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _charged;
}

void LEPFourJetsAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _charged;
}

void LEPFourJetsAnalysis::analyze(tEventPtr event, long, int, int ) {
  tPVector particles;

  if (_charged) {
    event->select(back_inserter(particles),
		  ThePEG::ParticleSelector<ChargedFinalState>());
  } else {
    event->select(back_inserter(particles),SelectFinalState());
  }


  //  copy fastjet particles from event record.  Templated fastjet
  //  method might leave units ambigouos.  Loop with integer index
  //  allows backtracing ThePEG particles if needed.
  vector<fastjet::PseudoJet> fastjet_particles;

  for (unsigned int j=0; j<particles.size(); j++) {
    fastjet::PseudoJet p(particles[j]->momentum().x()/GeV, 
			 particles[j]->momentum().y()/GeV, 
			 particles[j]->momentum().z()/GeV, 
			 particles[j]->momentum().e()/GeV);
    p.set_user_index(j);
    fastjet_particles.push_back(p);
  }
  
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm, 
				 recomb_scheme, strategy);
  fastjet::ClusterSequence cs(fastjet_particles, jet_def);

  vector<fastjet::PseudoJet> fastjets = cs.exclusive_jets_ycut(0.008);
  vector<fastjet::PseudoJet> sorted = fastjet::sorted_by_E(fastjets); 
  vector<Lorentz5Momentum> jets;

  if (sorted.size() == 4) {
    for (int j=0; j<4; ++j) {
      if ((cs.constituents(sorted[j])).size() == 1) {
	throw Exception() << "LEPFourJetsAnalysis: Trying to extract jet " 
			  << "momenta from a single particle." 
			  << Exception::warning;
      }
      LorentzMomentum newjet(sorted[j].px()*GeV, sorted[j].py()*GeV, 
			     sorted[j].pz()*GeV, sorted[j].e()*GeV);
      jets.push_back(newjet);
    }
    *_cchiBZ += abs(cosChiBZ(jets));
    *_cphiKSW += cosPhiKSW(jets);
    *_cthNR += abs(cosThetaNR(jets));
    *_ca34 += cosAlpha34(jets);
  }
}

ClassDescription<LEPFourJetsAnalysis> 
LEPFourJetsAnalysis::initLEPFourJetsAnalysis;
// Definition of the static class description member.

void LEPFourJetsAnalysis::Init() {

  static ClassDocumentation<LEPFourJetsAnalysis> documentation
    ("The LEP FourJets Analysis class",
     "The LEP FourJets analysis uses data from \\cite{Heister:2002tq}.",
     "%\\cite{Heister:2002tq}\n"
     "\\bibitem{Heister:2002tq}\n"
     "  A.~Heister {\\it et al.}  [ALEPH Collaboration],\n"
     "   ``Measurements of the strong coupling constant and the QCD colour factors\n"
     "  %using four-jet observables from hadronic Z decays,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 27}, 1 (2003).\n"
     "  %%CITATION = EPHJA,C27,1;%%\n"
     );

  static Switch<LEPFourJetsAnalysis,bool> interfaceChargedParticles
    ("ChargedParticles",
     "Wether or not to use charged particles only for this analysis",
     &LEPFourJetsAnalysis::_charged, true, false, false);
  static SwitchOption interfaceChargedParticlesYes
    (interfaceChargedParticles,
     "Yes",
     "Use charged particles only",
     true);
  static SwitchOption interfaceChargedParticlesNo
    (interfaceChargedParticles,
     "No",
     "Use all final state particles",
     false);

}

void LEPFourJetsAnalysis::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") 
    + name() + string(".top");
  ofstream output(fname.c_str());
  _ca34->normaliseToData();
  _cchiBZ->normaliseToData();
  _cphiKSW->normaliseToData();
  _cthNR->normaliseToData();
  // chisq
  double chisq,minfrac=0.05;
  unsigned int npoint;
  generator()->log() << "Output from LEPFourJetsAnalysis \n";
  _ca34->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for alpha_34 four jet distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  _cchiBZ->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for BZ four jet distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  _cphiKSW->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for KSW four jet distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  _cthNR->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for NR four jet distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  using namespace HistogramOptions;
  // output the plots
  _ca34->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "cosA0341 to DELPHI data",
		       "   GX  X               ",
		       "1/NdN/dcosA0341",
		       "          GX  X",
		       "cosA0341",
		       "   GX  X");
  _cchiBZ->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "cos|C0BZ1| to DELPHI data",
			 "    GX  X                ",
			 "1/NdN/dcos|C0BZ1|",
			 "           GX  X ",
			 "|cosC0BZ1|",
			 "    GX  X ");
  _cphiKSW->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "cosF0KSW1 to DELPHI data",
			 "   FX   X               ",
			 "1/NdN/dcosF0KSW1",
			 "          FX   X",
			 " cosF0KSW1",
			 "    FX   X");
  _cthNR->topdrawOutput(output,Frame|Errorbars|Ylog,
			 "RED",
			 "|cosQ0NR1| to DELPHI data",
			 "    GX  X                ",
			 "1/NdN/d|cosQ0NR1|",
			 "           GX  X ",
			 "|cosQ0NR1|",
			 "    GX  X ");
}

void LEPFourJetsAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // 4 jet angles 
  double ca34bins[] = {-1.0, -0.9, -0.8, -0.7, -0.6, 
		       -0.5, -0.4, -0.3, -0.2, -0.1, 
		       0.0, 0.1, 0.2, 0.3, 0.4, 
		       0.5, 0.6, 0.7, 0.8, 0.9, 
		       1.0};
  double ca34data[]={0.05279   ,0.04785   ,0.04541   ,0.04407   ,0.04461   ,
		     0.04225   ,0.04163   ,0.04166   ,0.04106   ,0.04180   ,
		     0.04171   ,0.04197   ,0.04193   ,0.04280   ,0.04245   ,
		     0.04412   ,0.03966   ,0.03199   ,0.01450   ,0.000417};
  double ca34errorstat[]={0.00079   ,0.00075   ,0.00073   ,0.00073   ,0.00074   ,
			  0.00072   ,0.00071   ,0.00071   ,0.00071   ,0.00072   ,
			  0.00071   ,0.00071   ,0.00071   ,0.00072   ,0.00071   ,
			  0.00074   ,0.00068   ,0.00058   ,0.00036   ,0.000062};
  double ca34errorsyst[]={0.00267 ,0.00236 ,0.00229 ,0.00227 ,0.00237 ,
			  0.00224 ,0.00214 ,0.00225 ,0.00215 ,0.00221 ,
			  0.00213 ,0.00212 ,0.00208 ,0.00207 ,0.00196 ,
			  0.00207 ,0.00163 ,0.00113 ,0.00033 ,0.000022};
  double ca34error[20];
  for(unsigned int ix=0;ix<20;++ix){ca34error[ix]=sqrt(sqr(ca34errorstat[ix])+
						       sqr(ca34errorsyst[ix]));}
  vector<double> bins  = vector<double>(ca34bins ,ca34bins +21);
  vector<double> data  = vector<double>(ca34data ,ca34data +20);
  vector<double> error = vector<double>(ca34error,ca34error+20);
  _ca34= new_ptr(Histogram(bins,data,error));
  double cchiBZbins[] = {0.00, 0.05, 0.10, 0.15, 0.20, 
			 0.25, 0.30, 0.35, 0.40, 0.45, 
			 0.50, 0.55, 0.60, 0.65, 0.70, 
			 0.75, 0.80, 0.85, 0.90, 0.95, 
			 1.00};
  double cchiBZdata[]={0.05455  ,0.05346  ,0.05450  ,0.05782  ,0.05730  ,
		       0.05880  ,0.05734  ,0.05973  ,0.06074  ,0.06387  ,
		       0.06483  ,0.06781  ,0.07144  ,0.07206  ,0.07887  ,
		       0.08601  ,0.09318  ,0.09865  ,0.11785  ,0.24115};
  double cchiBZerrorstat[]={0.00115  ,0.00113    ,0.00115  ,0.00120  ,0.00120  ,
			    0.00121  ,0.00118  ,0.00121  ,0.00122  ,0.00125  ,
			    0.00126  ,0.00128  ,0.00131  ,0.00129  ,0.00136  ,
			    0.00148  ,0.00146  ,0.00147  ,0.00159  ,0.00244};
  double cchiBZerrorsyst[]={0.00288,0.00287,0.00289,0.00306,0.00310,
			    0.00314,0.00299,0.00311,0.00319,0.00330,
			    0.00340,0.00332,0.00358,0.00329,0.00365,
			    0.00412,0.00406,0.00401,0.00432,0.01276};
  double cchiBZerror[20];
  for(unsigned int ix=0;ix<20;++ix){cchiBZerror[ix]=sqrt(sqr(cchiBZerrorstat[ix])+
							 sqr(cchiBZerrorsyst[ix]));}
  bins  = vector<double>(cchiBZbins ,cchiBZbins +21);
  data  = vector<double>(cchiBZdata ,cchiBZdata +20);
  error = vector<double>(cchiBZerror,cchiBZerror+20);
  _cchiBZ= new_ptr(Histogram(bins,data,error));
  double cphiKSWbins[] = {-1.0, -0.9, -0.8, -0.7, -0.6, 
			  -0.5, -0.4, -0.3, -0.2, -0.1, 
			  0.0, 0.1, 0.2, 0.3, 0.4, 
			  0.5, 0.6, 0.7, 0.8, 0.9, 
			  1.0};
  double cphiKSWdata[]={0.06378  ,0.03897  ,0.03558  ,0.03637  ,0.03597  ,
			0.03759  ,0.03696  ,0.03886  ,0.03801  ,0.03783  ,
			0.03342  ,0.03096  ,0.03033  ,0.02974  ,0.02976  ,
			0.02979  ,0.03068  ,0.03399  ,0.04234  ,0.09341};
  double cphiKSWerrorstat[]={0.00091  ,0.00066  ,0.00062  ,0.00065  ,0.00065  ,
			     0.00067  ,0.00065  ,0.00068  ,0.00065  ,0.00064  ,
			     0.00060  ,0.00059  ,0.00060  ,0.00060  ,0.00061  ,
			     0.00061  ,0.00062  ,0.00065  ,0.00072  ,0.00109};
  double cphiKSWerrorsyst[]={0.00362,0.00161,0.00143,0.00158,0.00157,
			     0.00172,0.00165,0.00177,0.00157,0.00143,
			     0.00130,0.00141,0.00154,0.00158,0.00172,
			     0.00172,0.00172,0.00191,0.00224,0.00560};
  double cphiKSWerror[20];
  for(unsigned int ix=0;ix<20;++ix){cphiKSWerror[ix]=sqrt(sqr(cphiKSWerrorstat[ix])+
							  sqr(cphiKSWerrorsyst[ix]));}
  bins  = vector<double>(cphiKSWbins ,cphiKSWbins +21);
  data  = vector<double>(cphiKSWdata ,cphiKSWdata +20);
  error = vector<double>(cphiKSWerror,cphiKSWerror+20);
  _cphiKSW= new_ptr(Histogram(bins,data,error));
  double cthNRbins[] = {0.00, 0.05, 0.10, 0.15, 0.20, 
			0.25, 0.30, 0.35, 0.40, 0.45, 
			0.50, 0.55, 0.60, 0.65, 0.70, 
			0.75, 0.80, 0.85, 0.90, 0.95, 
			1.00};
  double cthNRerror[20];
  double cthNRdata[]={0.06131  ,0.05888  ,0.05937  ,0.06104  ,0.05949  ,
		      0.06317  ,0.06632  ,0.06712  ,0.07040  ,0.07274  ,
		      0.07605  ,0.07707  ,0.08350  ,0.08779  ,0.08856  ,
		      0.09567  ,0.09632  ,0.10124  ,0.10139  ,0.12596};
  double cthNRerrorstat[]={0.00119  ,0.00114  ,0.00115  ,0.00118  ,0.00115  ,
			   0.00121  ,0.00125  ,0.00125  ,0.00129  ,0.00131  ,
			   0.00135  ,0.00135  ,0.00142  ,0.00147  ,0.00145  ,
			   0.00154  ,0.00154  ,0.00160  ,0.00158  ,0.00178};
  double cthNRerrorsyst[]={0.00281,0.00252,0.00272,0.00281,0.00271,
			   0.00292,0.00321,0.00324,0.00343,0.00354,
			   0.00377,0.00376,0.00422,0.00436,0.00431,
			   0.00503,0.00505,0.00534,0.00506,0.00630};
  for(unsigned int ix=0;ix<20;++ix){cthNRerror[ix]=sqrt(sqr(cthNRerrorstat[ix])+
							sqr(cthNRerrorsyst[ix]));}
  bins  = vector<double>(cthNRbins ,cthNRbins +21);
  data  = vector<double>(cthNRdata ,cthNRdata +20);
  error = vector<double>(cthNRerror,cthNRerror+20);
  _cthNR= new_ptr(Histogram(bins,data,error));
}

