// -*- C++ -*-
//
// ShowerInteraction.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerInteraction_H
#define HERWIG_ShowerInteraction_H

namespace Herwig {

/** \ingroup Shower
 *
 *  Handy header file to be included in all Shower classes.
 *  It contains only some useful enums.
 */

  namespace ShowerInteraction {

    /**
     *  Enum for the type of interaction
     */
    enum Type { UNDEFINED=-1, QCD, QED, QEDQCD, EW, ALL };
  }

  namespace ShowerPartnerType {
    /**
     *  Enum for the type of shower partner
     */
    enum Type {Undefined,QCDColourLine,QCDAntiColourLine,QED,EW};
  }

}
#endif // HERWIG_ShowerInteraction_H
