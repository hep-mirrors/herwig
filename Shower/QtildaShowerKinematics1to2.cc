// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/CLHEPWrap/Lorentz5Vector.h"
#include "Pythia7/Repository/CurrentGenerator.h"

using namespace Herwig;


QtildaShowerKinematics1to2::~QtildaShowerKinematics1to2() {}

vector<Lorentz5Momentum> QtildaShowerKinematics1to2::getBasis() {

  //cerr << "*";  
  vector<Lorentz5Momentum> dum; 
  //cerr << "*";
  dum.clear(); 
  //cerr << "*";
  dum.push_back( _pVector ); 
  //cerr << "*";
  dum.push_back( _nVector ); 
  //cerr << "*";
  //cerr << dum[0] << endl; 
  //cerr << dum[1] << endl; 
  return dum; 
}

// --- protected tools ---

Lorentz5Momentum QtildaShowerKinematics1to2::
sudakov2Momentum(double alpha, double beta, Energy px, Energy py) {

  // gives loads of output to check the transformations... 

  const Hep3Vector beta_bb = -(_pVector + _nVector).boostVector();

  // see also these methods!
  //    Hep3Vector findBoostToCM() const;
  //    // Boost needed to get to center-of-mass  frame:
  //            // w.findBoostToCM() == - w.boostVector()
  //            // w.boost(w.findBoostToCM()) == w.rest4Vector()
  
  //    Hep3Vector findBoostToCM( const HepLorentzVector & w ) const;
  //    // Boost needed to get to combined center-of-mass frame:
  //            // w1.findBoostToCM(w2) == w2.findBoostToCM(w1)
  //            // w.findBoostToCM(w) == w.findBoostToCM()

  Lorentz5Momentum p_bb = _pVector;
  Lorentz5Momentum n_bb = _nVector; 
  p_bb.boost( beta_bb );
  n_bb.boost( beta_bb );

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    if ( (p_bb[0] + n_bb[0] > 1.0*eV) || (p_bb[1] + n_bb[1] > 1.*eV) || (p_bb[2] + n_bb[2] > 1.*eV) ) {
      CurrentGenerator::log() << "QtildaShowerKinematics1to2::sudakov2Momentum(): " 
			      << "  Warning! check failed: haven't got to b2b-frame!"
			      << endl;  
    }
  }  

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "QtildaShowerKinematics1to2::sudakov2Momentum(): " 
			    << "==> start <==="
			    << endl;  
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "  called with (a, b, px, py) = (" << alpha
			    << ", " << beta
			    << ", " << px
			    << ", " << py
			    << ") " 
			    << endl
			    << "  LAB: p = " << _pVector 
			    << endl 
			    << "       n = " << _nVector << endl 
			    << "       p+n = " << _pVector + _nVector << endl 
			    << "       (p.n, p2, n2) = (" << _pVector*_nVector
			    << ", " << _pVector.mass2() 
			    << ", " << _nVector.mass2()
			    << ")" 
			    << endl
			    << "  to b2b frame: beta = " << beta_bb
			    << endl 
			    << "  b2b: p = " << p_bb 
			    << endl 
			    << "       n = " << n_bb 
			    << endl
			    << "       p+n = " << p_bb + n_bb
			    << endl 
			    << "       (p.n, p2, n2) = (" << p_bb*n_bb
			    << ", " << p_bb.mass2() 
			    << ", " << n_bb.mass2()
			    << ")"
			    << endl;
  }

  // set first in b2b frame along z-axis (assuming that p and n are
  // b2b as checked above)
  Lorentz5Momentum dq(0.0, 0.0, (alpha - beta)*p_bb.vect().mag(), 
     alpha*p_bb[3] + beta*n_bb[3] );

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "  q = " << dq << " ...(= a p'' + b n'') along z-axis" 
			    << endl
			    << "  zxs: in (theta, phi) = (" << p_bb.theta()
			    << ", " << p_bb.phi() << ")-direction. "
			    << endl
			    << "       p_bb/|p_bb| = " << p_bb.vect()/p_bb.vect().mag()
			    << endl
			    << "       q3/|q3|     = " << dq.vect().rotateUz( p_bb.vect()/p_bb.vect().mag() )/dq.vect().rotateUz( p_bb.vect()/p_bb.vect().mag() ).mag()
			    << endl; 
				 }

  // add transverse components
  dq.setPx(px);
  dq.setPy(py);

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "  ...including p_T:" 
			    << endl
			    << "  q = " << dq 
			    << endl 
			    << "  q2 = " << sqr(dq) 
			    << ", |q3| = " << dq.vect().mag()
			    << endl; 
  }
  
  // rotate to have z-axis parallel to p
  dq.rotateUz( p_bb.vect()/p_bb.vect().mag() );

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    Hep3Vector qpar = (((dq.vect())*(p_bb.vect()))/sqr((p_bb.vect()).mag()))*p_bb.vect();
    Hep3Vector qperp = dq.vect() - qpar;

    CurrentGenerator::log() << "  ...after rotation:"
			    << endl 
			    << "  q = " << dq 
			    << ", q2 = " << sqr(dq) 
			    << endl 
			    << "  q-drctn (theta, phi) = (" 
			    << qpar.theta()
			    << ", " 
			    << qpar.phi() 
			    << ")" 
			    << endl 
			    << "  |q3|/given = " 
			    << qpar.mag() 
			    << "/" << (alpha - beta)*p_bb.vect().mag()
			    << ", perp/given = " << qperp.mag() 
			    << "/" << sqrt(sqr(px) + sqr(py)) 
			    << endl;
  }

  // boost back 
  dq.boost( -beta_bb ); 
  dq.rescaleMass(); 

  // check consistency by getting back the Sudakov components from
  // the constructed momentum in the given basis 
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    double da = (dq*_nVector)/(_pVector*_nVector);  
    double db = (dq*_pVector - da*_pVector.m2())/(_pVector*_nVector);  
    double dq2 = - sqr( dq - da*_pVector - db*_nVector );

    CurrentGenerator::log() << "  ...after boost:" 
			    << "  q = " << dq 
			    << ", q2 = " << sqr(dq) 
			    << endl
			    << "  check \t given \t\t result \t ok?" << endl 
			    << "  alpha \t " 
			    << alpha << " \t " << da
			    << (da - alpha < 1e-10 ? " \t ok " : " \t not ok!") 
			    << endl 
			    << "  beta  \t " << beta << " \t " << db 
			    << (db - beta < 1e-10 ? " \t ok " : " \t not ok!")
			    << endl 
			    << "  perp2 \t " 
			    << dq2 << " \t " << sqr(px) + sqr(py) 
			    << (sqr(px) + sqr(py) > 0 ? 
				( (dq2 - sqr(px) - sqr(py))/sqr(px + sqr(py)) 
				  < 1e-10 ? " \t ok " : " not \t ok!") : " \t can't check.")
			    << endl
			    << "  (m2, (am)^2, 2ab p.n, (a m)^2+2ab p.n, q2) = " << endl 
			    << "  given  (" 
			    << _pVector.m2() << ", " 
			    << sqr(alpha)*_pVector.m2() << ", "
			    << 2*alpha*beta*(_pVector*_nVector) << ", "
			    << sqr(alpha)*_pVector.m2() + 2*alpha*beta*(_pVector*_nVector) 
			    << ", "
			    << sqr(alpha)*_pVector.m2() + 2*alpha*beta*(_pVector*_nVector) - sqr(px) - sqr(py) 
			    << ")" << endl 
			    << "  result (" 
			    << _pVector.m2() << ", " 
			    << sqr(da)*_pVector.m2() << ", "
			    << 2*da*db*(_pVector*_nVector) << ", "
			    << sqr(da)*_pVector.m2() + 2*da*db*(_pVector*_nVector) 
			    << ", "
			    << sqr(da)*_pVector.m2() + 2*da*db*(_pVector*_nVector) - dq2 
			    << ")"
			    << endl; 
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "QtildaShowerKinematics1to2::sudakov2Momentum(): " 
			    << "==> end <==="
			    << endl;  
  }

  return dq; 
  
}

