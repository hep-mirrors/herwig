// -*- C++ -*-
//
// Remapper.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_Remapper_H
#define Herwig_Remapper_H
//
// This is the declaration of the Remapper class.
//

#include <iostream>
#include <map>
#include "Herwig/Utilities/XML/Element.h"

namespace Herwig {

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Remapper adapts indivdual MC dimensions
 */
struct Remapper {

  std::map<double,double> weights;

  struct SelectorEntry {

    double lower;
    double upper;
    double value;

  };

  std::map<double,SelectorEntry> selector;

  double minSelection;

  bool smooth;

  Remapper();

  Remapper(unsigned int nBins,
	   double nMinSelection,
	   bool nSmooth);

  void fill(double x, double w);

  void finalize();

  std::pair<double,double> generate(double r) const;

  void fromXML(const XML::Element& elem);

  XML::Element toXML() const;

  void test(size_t n, std::ostream&);

};

}

#endif // Herwig_Remapper_H

