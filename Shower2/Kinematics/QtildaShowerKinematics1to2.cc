// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QtildaShowerKinematics1to2 class.
//

#include "QtildaShowerKinematics1to2.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "QtildaShowerKinematics1to2.tcc"
#endif


using namespace Herwig;

QtildaShowerKinematics1to2::~QtildaShowerKinematics1to2() {}

vector<Lorentz5Momentum> QtildaShowerKinematics1to2::getBasis() {
  vector<Lorentz5Momentum> dum;
  dum.push_back( _pVector );
  dum.push_back( _nVector );
  return dum; 
}

Lorentz5Momentum QtildaShowerKinematics1to2::
sudakov2Momentum(double alpha, double beta, Energy px, Energy py,unsigned int iopt) 
{
  Lorentz5Momentum dq;
  if(iopt==0)
    {
      const Hep3Vector beta_bb = -(_pVector + _nVector).boostVector();
      Lorentz5Momentum p_bb = _pVector;
      Lorentz5Momentum n_bb = _nVector; 
      p_bb.boost( beta_bb );
      n_bb.boost( beta_bb );
      // set first in b2b frame along z-axis (assuming that p and n are
      // b2b as checked above)
      dq=Lorentz5Momentum(0.0, 0.0, (alpha - beta)*p_bb.vect().mag(), 
			  alpha*p_bb.t() + beta*n_bb.t());
      // add transverse components
      dq.setPx(px);
      dq.setPy(py);
      // rotate to have z-axis parallel to p
      dq.rotateUz( p_bb.vect()/p_bb.vect().mag() );
      // boost back 
      dq.boost( -beta_bb ); 
      dq.rescaleMass(); 
      // return the momentum
    }
  else
    {
      const Hep3Vector beta_bb = -pVector().boostVector();
      Lorentz5Momentum p_bb = pVector();
      Lorentz5Momentum n_bb = nVector(); 
      p_bb.boost( beta_bb );
      n_bb.boost( beta_bb );
      // set first in b2b frame along z-axis (assuming that p and n are
      // b2b as checked above)
      dq=Lorentz5Momentum (0.0, 0.0, 0.5*beta*pVector().mass(), 
			  alpha*pVector().mass() + 0.5*beta*pVector().mass());
      // add transverse components
      dq.setPx(px);
      dq.setPy(py);
      // rotate to have z-axis parallel to n
      dq.rotateUz( n_bb.vect().unit());
      // boost back 
      dq.boost( -beta_bb ); 
      dq.rescaleMass();
    }
  return dq; 
}
