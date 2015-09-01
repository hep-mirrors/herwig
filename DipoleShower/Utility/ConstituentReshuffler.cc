// -*- C++ -*-
//
// ConstituentReshuffler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ConstituentReshuffler class.
//

#include <config.h>
#include "ConstituentReshuffler.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include <limits>

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "DipolePartonSplitter.h"

#include "Herwig/Utilities/GSLBisection.h"

#include "Herwig/DipoleShower/DipoleShowerHandler.h"

using namespace Herwig;

ConstituentReshuffler::ConstituentReshuffler() 
  : HandlerBase() {}

ConstituentReshuffler::~ConstituentReshuffler() {}

IBPtr ConstituentReshuffler::clone() const {
  return new_ptr(*this);
}

IBPtr ConstituentReshuffler::fullclone() const {
  return new_ptr(*this);
}

double ConstituentReshuffler::ReshuffleEquation::aUnit() {
  return 1.;
}

double ConstituentReshuffler::ReshuffleEquation::vUnit() {
  return 1.;
}

double ConstituentReshuffler::ReshuffleEquation::operator() (double xi) const {

  double r = - w/GeV;

  for (PList::iterator p = p_begin; p != p_end; ++p) {
    r += sqrt(sqr((**p).dataPtr()->constituentMass()) +
	      xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))) / GeV;
  }

  return r;  

}

void ConstituentReshuffler::reshuffle(PList& out,
				      PPair& in,
				      PList& intermediates) {

  if (out.size() == 0)
    return;

  if (out.size() == 1) {

    PPtr recoiler;
    PPtr parton = out.front();

    if (DipolePartonSplitter::colourConnected(parton,in.first) &&
	DipolePartonSplitter::colourConnected(parton,in.second)) {
      if (UseRandom::rnd() < .5)
	recoiler = in.first;
      else
	recoiler = in.second;
    } else if (DipolePartonSplitter::colourConnected(parton,in.first)) {
      recoiler = in.first;
    } else if (DipolePartonSplitter::colourConnected(parton,in.second)) {
      recoiler = in.second;
    } else assert(false);

    assert(abs(recoiler->momentum().vect().perp2()/GeV2) < 1e-6);

    double sign = recoiler->momentum().z() < 0.*GeV ? -1. : 1.;

    Energy2 qperp2 = parton->momentum().perp2();

    if (qperp2/GeV2 < Constants::epsilon) {
      // no emission off a 2 -> singlet process which
      // needed a single forced splitting: should never happen (?)
      assert(false);
      throw Veto();
    }

    Energy2 m2 = sqr(parton->dataPtr()->constituentMass());
      
    Energy abs_q = parton->momentum().vect().mag();
    Energy qz = parton->momentum().z();
    Energy abs_pz = recoiler->momentum().t();
    assert(abs_pz > 0.*GeV);

    Energy xi_pz = sign*(2.*qperp2*abs_pz + m2*(abs_q + sign*qz))/(2.*qperp2);
    Energy x_qz = (2.*qperp2*qz + m2*(qz+sign*abs_q))/(2.*qperp2);

    Lorentz5Momentum recoiler_momentum 
      (0.*GeV,0.*GeV,xi_pz,xi_pz < 0.*GeV ? - xi_pz : xi_pz);

    recoiler_momentum.rescaleMass();

    Lorentz5Momentum parton_momentum 
      (parton->momentum().x(),parton->momentum().y(),x_qz,sqrt(m2+qperp2+x_qz*x_qz));

    parton_momentum.rescaleMass();

    PPtr n_parton = new_ptr(Particle(parton->dataPtr()));
    n_parton->set5Momentum(parton_momentum);

    DipolePartonSplitter::change(parton,n_parton,false);

    out.pop_front();
    intermediates.push_back(parton);
    out.push_back(n_parton);

    PPtr n_recoiler = new_ptr(Particle(recoiler->dataPtr()));
    n_recoiler->set5Momentum(recoiler_momentum);

    DipolePartonSplitter::change(recoiler,n_recoiler,true);

    intermediates.push_back(recoiler);

    if (recoiler == in.first) {
      in.first = n_recoiler;
    }

    if (recoiler == in.second) {
      in.second = n_recoiler;
    }

    return;

  }

  Energy zero (0.*GeV);
  Lorentz5Momentum Q (zero,zero,zero,zero);
    
  for (PList::iterator p = out.begin();
       p != out.end(); ++p) {
    Q += (**p).momentum();
  }

  Boost beta = Q.findBoostToCM();

  list<Lorentz5Momentum> mbackup;

  bool need_boost = (beta.mag2() > Constants::epsilon);

  if (need_boost) {

    for (PList::iterator p = out.begin();
	 p != out.end(); ++p) {
      Lorentz5Momentum mom = (**p).momentum();
      mbackup.push_back(mom);
      (**p).set5Momentum(mom.boost(beta));
    }

  }

  double xi;

  ReshuffleEquation solve (Q.m(),out.begin(),out.end());

  GSLBisection solver(1e-10,1e-8,10000);

  try {
    xi = solver.value(solve,0.0,1.1);
  } catch (GSLBisection::GSLerror) {
    throw DipoleShowerHandler::RedoShower();
  } catch (GSLBisection::IntervalError) {
    throw DipoleShowerHandler::RedoShower();
  }

  PList reshuffled;

  list<Lorentz5Momentum>::const_iterator backup_it;
  if (need_boost)
    backup_it = mbackup.begin();

  for (PList::iterator p = out.begin();
       p != out.end(); ++p) {

    PPtr rp = new_ptr(Particle((**p).dataPtr()));

    DipolePartonSplitter::change(*p,rp,false);

    Lorentz5Momentum rm (xi*(**p).momentum().x(),
			 xi*(**p).momentum().y(),
			 xi*(**p).momentum().z(),
			 sqrt(sqr((**p).dataPtr()->constituentMass()) +
			      xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))));

    rm.rescaleMass();

    if (need_boost) {
      (**p).set5Momentum(*backup_it);
      ++backup_it;
      rm.boost(-beta);
    }

    rp->set5Momentum(rm);

    intermediates.push_back(*p);
    reshuffled.push_back(rp);
      

  }

  out.clear();
  out.splice(out.end(),reshuffled);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ConstituentReshuffler::persistentOutput(PersistentOStream &) const {
}

void ConstituentReshuffler::persistentInput(PersistentIStream &, int) {
}

ClassDescription<ConstituentReshuffler> ConstituentReshuffler::initConstituentReshuffler;
// Definition of the static class description member.

void ConstituentReshuffler::Init() {

  static ClassDocumentation<ConstituentReshuffler> documentation
    ("The ConstituentReshuffler class implements reshuffling "
     "of partons on their nominal mass shell to their constituent "
     "mass shells.");

}

