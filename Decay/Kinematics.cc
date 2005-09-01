// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Kinematics class.
//

#include "Kinematics.h"
#include <ThePEG/CLHEPWrap/Lorentz5Vector.h>
#include <ThePEG/CLHEPWrap/LorentzVector.h>
#include <ThePEG/CLHEPWrap/LorentzRotation.h>
#include <math.h>

using namespace Herwig;
using namespace ThePEG;


void MyKinematics::twoBodyDecay(const Lorentz5Momentum & p,  
			      const double m1, const double m2,
			      const Vector3 & unitDir1,
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) 
{
  if ( p.m() >= m1 + m2  &&  m1 >= 0.0  &&  m2 >= 0.0  ) {
    Momentum3 pstarVector = unitDir1;
    pstarVector *= CMMomentum(p.mass(),m1,m2);
    p1 = Lorentz5Momentum(m1,pstarVector);
    p2 = Lorentz5Momentum(m2,-pstarVector);
    p1.boost( p.boostVector() );   // boost from CM to LAB
    p2.boost( p.boostVector() );
  }
}

void MyKinematics::twoBodyDecay(const Lorentz5Momentum & p,                         
			      const double m1, const double m2,                 
			      const double cosThetaStar1, const double phiStar1, 
			      Lorentz5Momentum & p1, Lorentz5Momentum & p2 ) 
{ twoBodyDecay(p,m1,m2,unitDirection(cosThetaStar1,phiStar1),p1,p2); }

/* Doesn't work with new LorentzRotation and isn't used anyway PR 1/9/05
    
     * Create matrix that rotates a 5-momentum  to z-axis, then rotates by \f$\phi\f$ 
     * (where \f$\phi\f$ is given as two double args, \f$\cos\phi\f$ and \f$\sin\phi\f$.
     * @param mom The 5-momentum
     * @param cos Cosine of the angle  \f$\phi\f$, \f$\cos\phi\f$.
     * @param sin Sine of the angle  \f$\phi\f$, \f$\sin\phi\f$.
     
      static LorentzRotation rotation(const Lorentz5Momentum mom, const double cos,
  				    const double sin);
// Utility to create a rotation matrix to take p to the z axis, rotate it
// by psi about the z axis (cp = cos(psi), sp = sin(psi))
LorentzRotation MyKinematics::rotation(const Lorentz5Momentum p,
				     const double cp, 
				     const double sp) 
{
   double pt, pp, ct, st, cf, sf;

#define PTCUT 1E-20
#define sq(a) (a)*(a)

   pt = sq(p.x())+ sq(p.y());
   pp = sq(p.z()) + pt;
   if(pt < pp*PTCUT) {
      ct = p.z() > 0 ? 1.0 : -1.0;
      st = 0.0; cf = 1.0; sf = 0.0;
   } else {
      pp = sqrt(pp); pt = sqrt(pt); 
      ct = p.z()/pp; st = pt/pp; cf = p.x()/pt; sf = p.y()/pt;
   }

#undef PTCUT
#undef sq

   return LorentzRotation(LorentzVector(cp*cf*ct+sp*sf, -cp*sf+sp*cf*ct, cf*st, 0),
                          LorentzVector(cp*sf*ct-sp*cf,  cp*cf+sp*sf*ct, sf*st, 0),
                          LorentzVector(-cp*st,          -sp*st,         ct,    0),
                          LorentzVector(0,               0,              0,     1));
}
*/

void MyKinematics::generateAngles(tEGPtr gen, double &ct, double &az) {
  ct = gen->rnd()*2.0 - 1.0;
  az = gen->rnd()*2.0*M_PI;
}
