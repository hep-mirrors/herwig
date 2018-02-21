// -*- C++ -*-
//
// KinematicHelpers.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2018 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KinematicHelpers_H
#define HERWIG_KinematicHelpers_H

namespace QTildeKinematics {



inline Energy2 pT2_FSR(Energy2 qt2, double z, Energy2 m02, Energy2 m12, Energy2 m22) {
	return sqr(z*(1-z))*(qt2-m02) - m12*(1-z) - m22*z;
}


inline Energy2 pT2_ISR(Energy2 qt2, double z, Energy2 m2) {
	return sqr(1-z)*qt2 - m2*z;
}



inline Energy pT_FSR(Energy2 qt2, double z, Energy2 m02, Energy2 m12, Energy2 m22) {
	return sqrt( pT2_FSR(qt2,z,m02,m12,m22) );
}


inline Energy pT_ISR(Energy2 qt2, double z, Energy2 m2) {
	return sqrt( pT2_ISR(qt2,z,m2) );
}




}



#endif
