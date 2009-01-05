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

using namespace Herwig;

void y23::analyze(tEventPtr event, long, int, int ) {
  _kint.clearMap();
  tPVector particles;
  event->selectFinalState(back_inserter(particles));

  //luc method was added to feb../KtJet../Distance.*
  KtJet::KtDistance* distance_scheme = new KtJet::KtDistanceLuc(1);
  
  KtJet::KtEvent ev = KtJet::KtEvent(_kint.convert(particles), 1, 1, 1);

  KtJet::KtEvent ev_luc = KtJet::KtEvent(_kint.convert(particles), 1, distance_scheme, 1);

  double yval =  ev.getYMerge(2);
  double yval_luc =  ev_luc.getYMerge(2);

  // ynm distributions
  *_y23 += log10( yval ); 
  *_y23_luc += log10( yval_luc ); 
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
 
  string fname2 = generator()->filename() + string("-")
    + name() + string("luc.top");
  ofstream output2(fname2.c_str());

  using namespace HistogramOptions;

  _y23->normaliseToCrossSection();
  _y23_luc->normaliseToCrossSection();

  _y23->topdrawOutput(output,Frame|Errorbars,
		      "RED",
		      "y23",
		      "",
		      "",
		      "",
		       "log10(y23)",
		       "");

  _y23_luc->topdrawOutput(output2,Frame|Errorbars,
		      "RED",
		      "y23 luc",
		      "",
		      "",
		      "",
		      "log10(y23)",
		      "");
}

void y23::doinitrun() {
  AnalysisHandler::doinitrun();
  _y23 = new_ptr(Histogram(-5., 0., 100 ));
  _y23_luc = new_ptr(Histogram(-5., 0., 100 ));
}

using namespace KtJet;
//luc functions

KtDistanceLuc::KtDistanceLuc(int collision_type) : m_type(collision_type), m_name("luclus") {}
  //KtDistanceAngle::~KtDistanceAngle() {}
std::string KtDistanceLuc::name() const {return m_name;}

KtFloat KtDistanceLuc::operator()(const KtLorentzVector & a) const {
  KtFloat kt, r, costh;
  const KtFloat small = 0.0001;     // ??? Should be defined somewhere else?
  switch (m_type) {            // direction of beam depends on collision type
  case 1:
    return -1;               // e+e- : no beam remnant, so result will be ignored anyway
    break;
  case 2:                    // ep (p beam -z direction)
    costh = -(a.cosTheta());
    break;
  case 3:                    // pe (p beam +z direction)
    costh = a.cosTheta();
    break;
  case 4:                    // pp (p beams in both directions)
    costh = fabs(a.cosTheta());
    break;
  default:                   // type out of range - WARNING ???
    costh = 0.;
    break;
  }
  r = 2*(1-costh);
  if (r<small) r = a.perp2()/a.vect().mag2();  // Use approx if close to beam
  kt = a.e()*a.e() * r;
  return kt;
}

KtFloat KtDistanceLuc::operator()(const KtLorentzVector & a, const KtLorentzVector & b) const {
  KtFloat esq = a.e() * b.e() * a.e() * b.e() 
    /( a.e() + b.e() ) / ( a.e() + b.e() );
  KtFloat costh = a.vect().cosTheta(b.vect());
  return 2. * esq * (1 - costh);
}


