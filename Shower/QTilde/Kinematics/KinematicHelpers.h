// -*- C++ -*-
//
// KinematicHelpers.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2018-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KinematicHelpers_H
#define HERWIG_KinematicHelpers_H

namespace QTildeKinematics {



inline Energy2 pT2_FSR(Energy2 qt2, double z, Energy2 m02, Energy2 m12, Energy2 m22,
		       Energy2 q12, Energy2 q22) {
  const double z1z = z*(1.-z);
  return z1z*(z1z*qt2 + m02 -m12-m22) - q12*sqr(1.-z) - q22*sqr(z);
  //return z1z*(z1z*qt2 + m02) - q12*(1.-z) - q22*z;
}

inline Energy2 pT2_ISR(Energy2 qt2, double z, Energy2 m22) {
  return sqr(1.-z)*qt2 - m22*z;
}

inline Energy2 pT2_Decay(Energy2 qt2, double z, Energy2 m02, Energy2 m22) {
  return sqr(1.-z)*(qt2 - m02) - m22*z;
}



inline Energy pT_FSR(Energy2 qt2, double z, Energy2 m02, Energy2 m12, Energy2 m22,
		       Energy2 q12, Energy2 q22) {
  return sqrt( pT2_FSR(qt2,z,m02,m12,m22,q12,q22) );
}


inline Energy pT_ISR(Energy2 qt2, double z, Energy2 m22) {
  return sqrt( pT2_ISR(qt2,z,m22) );
}

inline Energy pT_Decay(Energy2 qt2, double z, Energy2 m02, Energy2 m22) {
  return sqrt( pT2_Decay(qt2,z,m02,m22) );
}


inline Energy2 q2_FSR(Energy2 pt2, double z, Energy2 m12, Energy2 m22) {
  return m12/z + m22/(1.-z) + pt2/z/(1.-z);
}

// inline Energy2 q2_ISR(Energy2 pt2, double z, Energy2 m22) {
//   return m22/(1.-z) + pt2/z/(1.-z);
// }


}

#endif
