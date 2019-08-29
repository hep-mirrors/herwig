// -*- C++ -*-
//
// IsoSpin.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IsoSpin_H
#define HERWIG_IsoSpin_H
//
// This is the declaration of the Isospin namespace and values.
//
namespace Herwig {
  namespace IsoSpin {
    
    /**
     *   Enum for the total isospin of the system
     */
    enum IsoSpin { IUnknown, IZero, IHalf, IOne};
    
    /**
     *   Third component
     */
    enum I3     { I3Unknown, I3MinusOne, I3MinusHalf, I3Zero, I3Half, I3One};
    
  }
}

#endif /* HERWIG_IsoSpin_H */
