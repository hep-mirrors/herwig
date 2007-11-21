// -*- C++ -*-
//
// Herwig.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_H
#define HERWIG_H
//
// This is the declaration of the <!id>Herwig.h<!!id> header file.
//
// CLASSDOC SUBSECTION Description:
//
// This header file should be included in all Herwig++ classes. <BR>
// At the moment, it is used mainly to define the macro HERWIG_NO_DEBUG <BR> 
// which switches off all the debugging information: see below <BR>
//      <I>  #define HERWIG_NO_DEBUG  YES  </I><BR>
// (alternatively, one could define such macro by compiling with the flag: <BR> 
//        g++ -c -DHERWIG_NO_DEBUG ) <BR>
//

#include "ThePEG/Config/ThePEG.h"

//***LOOKHERE***  #define HERWIG_NO_DEBUG  1

// Debugging in Herwig may be switched off completely by this compilation 
// switched, eliminating possible overhead in error checking.
#ifndef HERWIG_NO_DEBUG
#define HERWIG_DEBUG_LEVEL Herwig::HwDebug::level
#else
#define HERWIG_DEBUG_LEVEL 0
#endif

#endif /* HERWIG_H */

