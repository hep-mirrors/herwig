// -*- C++ -*-
//
// Reshuffler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Reshuffler_H
#define HERWIG_Reshuffler_H
//
// This is the declaration of the Reshuffler class.
//

#include "ThePEG/Config/ThePEG.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \author Simon Platzer, Stephen Webster
 * 
 * \brief The Reshuffler class implements reshuffling
 * of partons on their nominal mass shell to their constituent 
 * mass shells.
 *
 */
class Reshuffler {

protected:

  /**
   * The function object defining the equation
   * to be solved.
   */
  template<class PIterator, class MIterator>
  struct ReshuffleEquation {

    ReshuffleEquation (Energy q,
		       PIterator pp_begin,
		       PIterator pp_end,
		       MIterator mm_begin,
		       MIterator mm_end)
      : w(q), p_begin(pp_begin), p_end(pp_end),
	m_begin(mm_begin), m_end(mm_end) {}

    typedef double ArgType;
    typedef double ValType;

    static double aUnit() { return 1.; }
    static double vUnit() { return 1.; }

    double operator() (double xi) const {
      double r = - w/GeV;
      PIterator p = p_begin;
      MIterator m = m_begin;
      for (; p != p_end; ++p, ++m) {
	r += sqrt(sqr(*m) +
		  xi*xi*(sqr((**p).momentum().t())-sqr((**p).dataPtr()->mass()))) / GeV;
      }
      return r;  
    }

    Energy w;

    PIterator p_begin;
    PIterator p_end;

    MIterator m_begin;
    MIterator m_end;

  };

  /**
   * Reshuffle to consitutent masses
   */
  void reshuffle(const PVector& particles,
		 const vector<Energy>& masses) const;

};

}

#endif /* HERWIG_Reshuffler_H */

