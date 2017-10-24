// -*- C++ -*-
//
// ShowerInteraction.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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

/**
 *  Enum for the type of interaction
 */
enum class ShowerInteraction { 
	UNDEFINED=-1, 
	QCD, 
	QED, 
	Both 
};

/**
 *  Enum for the type of shower partner
 */
enum class ShowerPartnerType {
	Undefined,
	QCDColourLine,
	QCDAntiColourLine,
	QED
};

}
#endif // HERWIG_ShowerInteraction_H
