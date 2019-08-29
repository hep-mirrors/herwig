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
  namespace Strangeness {
    /**
     *  Strange content
     */
    enum Strange { Unknown, ssbar, Zero, PlusOne, MinusOne};
  }
  namespace Charm {
    /**
     *  Charm content
     */
    enum Charm { Unknown, ccbar, Zero, PlusOne, MinusOne};
  }
  namespace Beauty {
    /**
     *  Bottom content
     */
    enum Bottom { Unknown, bbbar, Zero, PlusOne, MinusOne};
  }

  struct FlavourInfo {

    /**
     *  Constructor
     */
    FlavourInfo() : I(IsoSpin::IUnknown), I3(IsoSpin::I3Unknown),
		    strange(Strangeness::Unknown), charm(Charm::Unknown), bottom(Beauty::Unknown)
    {}
    
    /**
     *  Constructor
     */
    FlavourInfo(IsoSpin::IsoSpin Iin, IsoSpin::I3 I3in, Strangeness::Strange Sin=Strangeness::Unknown,
		Charm::Charm Cin=Charm::Unknown, Beauty::Bottom Bin=Beauty::Unknown) :
      I(Iin), I3(I3in),strange(Sin), charm(Cin), bottom(Bin)
    {}

    /**
     *  Total isospin
     */
    IsoSpin::IsoSpin I;
    
    /**
     *  \f$I_3\f$
     */
    IsoSpin::I3 I3;

    /**
     * Strange
     */
    Strangeness::Strange strange;

    /**
     * Charm
     */
    Charm::Charm charm;

    /**
     * Strange
     */
    Beauty::Bottom bottom;
  };
}

#endif /* HERWIG_IsoSpin_H */
