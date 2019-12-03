// -*- C++ -*-
//
// DiagramDrawer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_DiagramDrawer_H
#define Herwig_DiagramDrawer_H

#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief DiagramDrawer draws ASCII output from Tree2toNDiagram
 * objects for diagnostic purposes.
 *
 */
struct DiagramDrawer {

  /**
   * Draw a diagram
   */
  static void drawDiag(ostream&,const Tree2toNDiagram&);  

};

}

#endif // Herwig_DiagramDrawer_H
