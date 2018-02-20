// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerBasis class.
//

#include "ShowerBasis.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;


// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<ShowerBasis,Base>
describeHerwigShowerBasis("Herwig::ShowerBasis", "libHerwig.so");

void ShowerBasis::setBasis(const Lorentz5Momentum & p,
			   const Lorentz5Momentum & n,
			   Frame inframe) {
  pVector_ = p;
  nVector_ = n;
  frame_ = inframe;
  Boost beta_bb;
  if(frame()==BackToBack) {
    beta_bb = -(pVector_ + nVector_).boostVector();
  }
  else if(frame()==Rest) {
    beta_bb = -pVector().boostVector();
  }
  else
    assert(false);
  Lorentz5Momentum p_bb = pVector();
  Lorentz5Momentum n_bb = nVector(); 
  p_bb.boost( beta_bb );
  n_bb.boost( beta_bb );
  // rotate to have z-axis parallel to p/n
  Axis axis;
  if(frame()==BackToBack) {
    axis = p_bb.vect().unit();
  }
  else if(frame()==Rest) {
    axis = n_bb.vect().unit();
  }
  else
    assert(false);
  LorentzRotation rot;
  if(axis.perp2()>1e-10) {
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    rot.rotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  }
  else if(axis.z()<0.) {
    rot.rotate(Constants::pi,Axis(1.,0.,0.));
  }
  xPerp_=LorentzVector<double>(1.,0.,0.,0.);
  yPerp_=LorentzVector<double>(0.,1.,0.,0.);
  xPerp_.transform(rot);
  yPerp_.transform(rot);
  // boost back 
  xPerp_.boost( -beta_bb );
  yPerp_.boost( -beta_bb );
}

void ShowerBasis::transform(const LorentzRotation & r) {
  pVector_ *= r;
  nVector_ *= r;
  xPerp_   *= r;
  yPerp_   *= r;
}

vector<Lorentz5Momentum> ShowerBasis::getBasis() const {
  vector<Lorentz5Momentum> dum;
  dum.push_back( pVector_ );
  dum.push_back( nVector_ );
  return dum; 
}
