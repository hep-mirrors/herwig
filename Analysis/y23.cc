// -*- C++ -*-
//
// y23.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the y23 class.
//

#include "y23.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Utilities/HerwigStrategy.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

void y23::analyze(tEventPtr event, long, int, int ) {

  tPVector particles;
  event->selectFinalState(back_inserter(particles));
  
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

  double yval =  cs.exclusive_ymerge(2);
  // ynm distributions
  *_y23 += log10( yval ); 
}

NoPIOClassDescription<y23> y23::inity23;
// Definition of the static class description member.

void y23::Init() {
  static ClassDocumentation<y23> documentation
    ("There is no documentation for the y23 class");
}

void y23::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") 
    + name() + string(".top");
  ofstream output(fname.c_str());
 
  using namespace HistogramOptions;

  _y23->normaliseToCrossSection();
 
  _y23->topdrawOutput(output,Frame|Errorbars,
		      "RED",
		      "y23",
		      "",
		      "",
		      "",
		       "log10(y23)",
		       "");
}

void y23::doinitrun() {
  AnalysisHandler::doinitrun();
  _y23 = new_ptr(Histogram(-5., 0., 100 ));
}


