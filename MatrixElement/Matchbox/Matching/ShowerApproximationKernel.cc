// -*- C++ -*-
//
// ShowerApproximationKernel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximationKernel class.
//

#include <config.h>
#include "ShowerApproximationKernel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ShowerApproximationGenerator.h"

using namespace Herwig;

ShowerApproximationKernel::ShowerApproximationKernel() 
  : thePresampling(false), thePresamplingPoints(10000), 
    theMaxTry(100000), theFreezeGrid(500000),
    sampler(0), theDoCompensate(false) {}

ShowerApproximationKernel::~ShowerApproximationKernel() {}

IBPtr ShowerApproximationKernel::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerApproximationKernel::fullclone() const {
  return new_ptr(*this);
}

void ShowerApproximationKernel::showerApproximationGenerator(Ptr<ShowerApproximationGenerator>::tptr gen) {
  theShowerApproximationGenerator = gen;
}

Ptr<ShowerApproximationGenerator>::tptr ShowerApproximationKernel::showerApproximationGenerator() const {
  return theShowerApproximationGenerator;
}

const vector<bool>& ShowerApproximationKernel::sampleFlags() {
  if ( !theFlags.empty() )
    return theFlags;
  theFlags.resize(nDim(),false);
  for ( int k = nDimBorn(); 
	k < nDimBorn() + dipole()->nDimRadiation(); ++k )
    theFlags[k] = true;
  return theFlags;
}

const pair<vector<double>,vector<double> >& ShowerApproximationKernel::support() {
  if ( !theSupport.first.empty() )
    return theSupport;
  vector<double> l(nDim(),0.0);
  vector<double> u(nDim(),1.0);
  theSupport.first = l;
  theSupport.second = u;
  return theSupport;
}

const vector<double>& ShowerApproximationKernel::parameterPoint() {
  theLastParameterPoint.resize(nDim());
  copy(bornCXComb()->lastRandomNumbers().begin(),
       bornCXComb()->lastRandomNumbers().end(),
       theLastParameterPoint.begin());
  theLastParameterPoint[evolutionVariable()] = 1.;
  return theLastParameterPoint;
}

void ShowerApproximationKernel::startPresampling() { 
  thePresampling = true;
}

void ShowerApproximationKernel::stopPresampling() {
  showerApproximationGenerator()->restore(); 
  thePresampling = false;
}

double ShowerApproximationKernel::evaluate(const vector<double>& r) {

  if ( presampling() ) {

    theLastBornPoint.resize(nDimBorn());
    copy(r.begin(),r.begin()+nDimBorn(),theLastBornPoint.begin());
    if ( !showerApproximationGenerator()->generate(theLastBornPoint) )
      return 0.;

  }

  assert(dipole()->splitting());
  realXComb()->clean();
  dipole()->setXComb(realXComb());
  for ( vector<StdXCombPtr>::const_iterator t = tildeXCombs().begin();
	t != tildeXCombs().end(); ++t ) {
    (**t).clean();
    (**t).matrixElement()->setXComb(*t);
  }
  if ( !dipole()->generateKinematics(&r[nDimBorn()]) )
    return 0.;

  double jac = 
    showerApproximation()->showerInvertedTildeKinematics() ?
    showerApproximation()->showerInvertedTildeKinematics()->jacobian() :
    dipole()->invertedTildeKinematics()->jacobian();

  showerApproximation()->setBornXComb(bornXComb());
  showerApproximation()->setRealXComb(realXComb());
  showerApproximation()->setTildeXCombs(tildeXCombs());
  showerApproximation()->setDipole(dipole());
  showerApproximation()->checkCutoff();
  showerApproximation()->getShowerVariables();

  if ( !dipole()->isInShowerPhasespace() )
    return 0.;

  return showerApproximation()->me2() * jac;

}

double ShowerApproximationKernel::generate() {

  if ( !sampler ) {

    sampler = new ExponentialGenerator();
    sampler->sampling_parameters().maxtry = maxtry();
    sampler->sampling_parameters().presampling_points = presamplingPoints();
    sampler->sampling_parameters().freeze_grid = freezeGrid();
    sampler->docompensate(theDoCompensate);
    sampler->function(this);
    sampler->initialize();

  }

  double res = 0.;

  while (true) {
    try {
      res = sampler->generate();
    } catch (exsample::exponential_regenerate&) {
      continue;
    } catch (exsample::hit_and_miss_maxtry&) {
      throw MaxTryException();
    } catch (exsample::selection_maxtry&) {
      throw MaxTryException();
    } 
    break;
  }

  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ShowerApproximationKernel::persistentOutput(PersistentOStream & os) const {
  os << theDipole << theShowerApproximation << theBornXComb 
     << theRealXComb << theTildeXCombs << thePresampling 
     << thePresamplingPoints << theMaxTry << theFreezeGrid << theFlags 
     << theSupport << theShowerApproximationGenerator 
     << theLastParameterPoint << theLastBornPoint
     << (sampler ? true : false) << theDoCompensate;
  if ( sampler )
    sampler->put(os);
}

void ShowerApproximationKernel::persistentInput(PersistentIStream & is, int) {
  bool haveSampler;
  is >> theDipole >> theShowerApproximation >> theBornXComb 
     >> theRealXComb >> theTildeXCombs >> thePresampling 
     >> thePresamplingPoints >> theMaxTry >> theFreezeGrid >> theFlags 
     >> theSupport >> theShowerApproximationGenerator 
     >> theLastParameterPoint >> theLastBornPoint
     >> haveSampler >> theDoCompensate; 
  if ( haveSampler ) {
    sampler = new ExponentialGenerator();
    sampler->get(is);
    sampler->function(this);
    sampler->initialize();
  }
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ShowerApproximationKernel,HandlerBase>
  describeHerwigShowerApproximationKernel("Herwig::ShowerApproximationKernel", "Herwig.so");

void ShowerApproximationKernel::Init() {

  static ClassDocumentation<ShowerApproximationKernel> documentation
    ("ShowerApproximationKernel generates emissions according to a "
     "shower approximation entering a NLO matching.");

}

