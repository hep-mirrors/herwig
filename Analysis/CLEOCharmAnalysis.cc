// -*- C++ -*-
//
// CLEOCharmAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEOCharmAnalysis class.
//

#include "CLEOCharmAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void CLEOCharmAnalysis::analyze(tEventPtr event, long, int, int) {
  _s = (event->incoming().first ->momentum()+
	event->incoming().second->momentum()).m2();
  double weight = event->weight();
  set<tPPtr> particles;
  StepVector steps = event->primaryCollision()->steps();
  for ( StepVector::const_iterator it = steps.begin()+2;
	it != steps.end(); ++it ) {
    (**it).select(inserter(particles), ThePEG::AllSelector());
  }
  tPVector output;
  for(set<tPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long id = abs((*it)->id());
    if(id==ParticleID::Dplus   || id==ParticleID::Dstarplus    ||
       id==ParticleID::D0      || id==ParticleID::Dstar0       )
      output.push_back(*it);
  }
  analyze(output,weight);
}

void CLEOCharmAnalysis::analyze(tPPtr particle, double weight) {
  // Calls analyze() for each particle.
  double xp = particle->momentum().vect().mag()/
    sqrt(0.25*_s-sqr(particle->mass()));
  int id = abs(particle->id());
  if(id==ParticleID::Dstarplus) {
    _histDstarplus->addWeighted(xp,weight);
  }
  else if(id==ParticleID::Dstar0) {
    _histDstar0   ->addWeighted(xp,weight);
  }
  else if(id==ParticleID::D0) {
    _histD0       ->addWeighted(xp,weight);
  }
  else if(id==ParticleID::Dplus) {
    _histDplus    ->addWeighted(xp,weight);
  }
}

NoPIOClassDescription<CLEOCharmAnalysis> CLEOCharmAnalysis::initCLEOCharmAnalysis;
// Definition of the static class description member.

void CLEOCharmAnalysis::Init() {

  static ClassDocumentation<CLEOCharmAnalysis> documentation
    ("CLEO Charm meson analysis class",
     "The CLEO Charm meson analysis uses data from \\cite{Artuso:2004pj}.",
     "%\\cite{Artuso:2004pj}\n"
     "\\bibitem{Artuso:2004pj}\n"
     "  M.~Artuso {\\it et al.}  [CLEO Collaboration],\n"
     "  %``Charm meson spectra in $e^{+} e^{-}$ annihilation at 10.5-GeV c.m.e,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 70}, 112001 (2004)\n"
     "  [arXiv:hep-ex/0402040].\n"
     "  %%CITATION = PHRVA,D70,112001;%%\n"
     );

}

void CLEOCharmAnalysis::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  double chisq=0.,minfrac=0.05;
  unsigned int ndegrees;
  _histDstarplus->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for CLEO D*+ distribution\n";
  _histDstar0   ->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for CLEO D*0 distribution\n";
  _histD0       ->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for CLEO D0 distribution\n";
  _histDplus    ->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for CLEO D+ distribution\n";
  _histDstarplus->topdrawOutput(output,Frame|Errorbars,
				"RED",
				"D2*+3",
				" X  X",
				"1/SdS/dx0p1",
				"  G G   X X",
				"x0p1",
				" X X");
  _histDstar0->topdrawOutput(output,Frame|Errorbars,
			     "RED",
			     "D2*03",
			     " X  X",
			     "1/SdS/dx0p1",
			     "  G G   X X",
			     "x0p1",
			     " X X");
  _histD0->topdrawOutput(output,Frame|Errorbars,
			 "RED",
			 "D203",
			 " X X",
			 "1/SdS/dx0p1",
			 "  G G   X X",
			 "x0p1",
			 " X X");
  _histDplus->topdrawOutput(output,Frame|Errorbars,
			    "RED",
			    "D2+3",
			    " X X",
			    "1/SdS/dx0p1",
			    "  G G   X X",
			    "x0p1",
			    " X X");
}

void CLEOCharmAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  double vals[]={0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,
		 0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
  // data and errors for D+
  double dplusdata     []={ 161, 320, 356, 413, 693,
			    909,1042,1271,1357,1370,
			   1291,1129, 952, 694, 449, 223,  74};
  double dpluserror    []={  83,  92,  92,  94,  60,
			     56,  53,  38,  38,  36,
			     34,  29,  25,  19,  13,   7,   4};
  // data and errors for D0
  double d0data        []={ 173, 431, 529, 882,1156,
			   1670,2349,2822,3194,3475,
			   3371,3007,2549,2008,1383, 829, 339,  90};
  double d0error       []={ 109, 186, 209,  84,  83,
			     94, 110, 122,  56,  58,
			     56,  51,  46,  39,  30,  21,  11,   5};
  // data and errors for D*+
  double dstarplusdata []={ 146, 253, 348, 494, 624,
			    920,1108,1244,1286,1248,
			   1113, 932, 723, 531, 310, 119};
  double dstarpluserror[]={  86,  60,  60,  60,  46,
			     50,  32,  33,  32,  31,
			     29,  25,  21,  17,  12,   8};
  // data and errors for D*0
  double dstar0data    []={ 108, 292, 387, 425, 594, 
			    546, 897,1085,1162,1230,
			   1198,1055, 865, 694, 471, 289, 121};
  double dstar0error   []={ 121, 115, 111, 103,  98,
			     92, 108,  70,  65,  64,
			     60,  52,  45,  36,  27,  20,  15};
  // D+ histogram
  vector<double> bins,data,error;
  for(unsigned int ix=0;ix<4;++ix) {
    if(ix==0) {
      bins  = vector<double>(vals+1,vals+19);
      data  = vector<double>(dplusdata ,dplusdata +17);
      error = vector<double>(dpluserror,dpluserror+17);
    }
    else if(ix==1) {
      bins  = vector<double>(vals,vals+19);
      data  = vector<double>(d0data ,d0data +18);
      error = vector<double>(d0error,d0error+18);
    }
    else if(ix==2) {
      bins  = vector<double>(vals+2,vals+19);
      data  = vector<double>(dstarplusdata ,dstarplusdata +16);
      error = vector<double>(dstarpluserror,dstarpluserror+16);
    }
    else if(ix==3) {
      bins  = vector<double>(vals+1,vals+19);
      data  = vector<double>(dstar0data ,dstar0data +17);
      error = vector<double>(dstar0error,dstar0error+17);
    }
    double norm=0.;
    for(unsigned int iy=0;iy<data.size();++iy) norm+=data[iy];
    norm *= 0.05; 
    for(unsigned int iy=0;iy<data.size();++iy) {
      data [iy] /= norm;
      error[iy] /= norm;
    }
    if(ix==0)      _histDplus=new_ptr(Histogram(bins,data,error));
    else if(ix==1) _histD0=new_ptr(Histogram(bins,data,error));
    else if(ix==2) _histDstarplus=new_ptr(Histogram(bins,data,error));
    else if(ix==3) _histDstar0=new_ptr(Histogram(bins,data,error));
  }
}
