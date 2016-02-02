// -*- C++ -*-
//
// ShowerApproximationKernel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ShowerApproximationKernel_H
#define Herwig_ShowerApproximationKernel_H
//
// This is the declaration of the ShowerApproximationKernel class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/Sampling/exsample/exponential_generator.h"

namespace Herwig {

using namespace ThePEG;

class ShowerApproximationGenerator;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ShowerApproximationKernel generates emissions according to a
 * shower approximation entering a NLO matching.
 *
 */
class ShowerApproximationKernel: public HandlerBase {

public:

  /**
   * Exception to communicate sampler maxtry events.
   */
  struct MaxTryException {};

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ShowerApproximationKernel();

  /**
   * The destructor.
   */
  virtual ~ShowerApproximationKernel();
  //@}

public:

  /**
   * Set the XComb object describing the Born process
   */
  void setBornXComb(tStdXCombPtr xc) { theBornXComb = xc; }

  /**
   * Return the XComb object describing the Born process
   */
  tcStdXCombPtr bornCXComb() const { return theBornXComb; }

  /**
   * Return the XComb object describing the Born process
   */
  tStdXCombPtr bornXComb() const { return theBornXComb; }

  /**
   * Set the XComb object describing the real emission process
   */
  void setRealXComb(tStdXCombPtr xc) { theRealXComb = xc; }

  /**
   * Return the XComb object describing the real emission process
   */
  tcStdXCombPtr realCXComb() const { return theRealXComb; }

  /**
   * Return the XComb object describing the real emission process
   */
  tStdXCombPtr realXComb() const { return theRealXComb; }

  /**
   * Set the tilde xcomb objects associated to the real xcomb
   */
  void setTildeXCombs(const vector<StdXCombPtr>& xc) { theTildeXCombs = xc; }

  /**
   * Return the tilde xcomb objects associated to the real xcomb
   */
  const vector<StdXCombPtr>& tildeXCombs() const { return theTildeXCombs; }

  /**
   * Set the dipole in charge for the emission
   */
  void setDipole(Ptr<SubtractionDipole>::tptr dip) { theDipole = dip; }

  /**
   * Return the dipole in charge for the emission
   */
  Ptr<SubtractionDipole>::tptr dipole() const { return theDipole; }

  /**
   * Set the shower approximation.
   */
  void showerApproximation(Ptr<ShowerApproximation>::tptr app) { theShowerApproximation = app; }

  /**
   * Return the shower approximation.
   */
  Ptr<ShowerApproximation>::tptr showerApproximation() const { return theShowerApproximation; }

  /**
   * Set the shower approximation generator.
   */
  void showerApproximationGenerator(Ptr<ShowerApproximationGenerator>::tptr);

  /**
   * Return the shower approximation generator.
   */
  Ptr<ShowerApproximationGenerator>::tptr showerApproximationGenerator() const;

  /**
   * Generate the next emission
   */
  double generate();

public:

  /**
   * Set a pt cut on the dipole to generate the radiation
   */
  void ptCut(Energy pt) { dipole()->ptCut(pt); }

  /**
   * Return the number of random numbers
   * needed to sample this kernel.
   */
  int nDim() const {
    return 
      nDimBorn() +
      dipole()->nDimRadiation();
  }

  /**
   * Return the number of random numbers
   * needed to sample the Born process.
   */
  int nDimBorn() const {
    return bornCXComb()->lastRandomNumbers().size();
  }

  /**
   * Flag, which variables are free variables.
   */
  const vector<bool>& sampleFlags();

  /**
   * Return the support of the splitting kernel.
   * The lower bound on the first variable is
   * assumed to correspond to the cutoff on the
   * evolution variable.
   */
  const pair<vector<double>,vector<double> >& support();

  /**
   * Return the parameter point associated to the splitting
   * previously supplied through fixParameters.
   */
  const vector<double>& parameterPoint();

  /**
   * Indicate that presampling of this kernel
   * will be performed in the next calls to
   * evaluate until stopPresampling() is called.
   */
  void startPresampling();

  /**
   * Indicate that presampling of this kernel
   * is done until startPresampling() is called.
   */
  void stopPresampling();

  /**
   * Indicate that a veto with the given kernel value and overestimate has occured.
   */
  void veto(const vector<double>&, double, double) {
    /** use born and real xcombs in here to figure out what we need to reweight;
	it should have its kinematic variables completed at this step */
  }

  /**
   * Indicate that an accept with the given kernel value and overestimate has occured.
   */
  void accept(const vector<double>&, double, double) {
    /** use born and real xcombs in here to figure out what we need to reweight;
	it should have its kinematic variables completed at this step */
  }

  /**
   * Return true, if currently being presampled
   */
  bool presampling() const { return thePresampling; }

  /**
   * Return the number of points to presample this
   * splitting generator.
   */
  unsigned long presamplingPoints() const { return thePresamplingPoints; }

  /**
   * Return the maximum number of trials
   * to generate a splitting.
   */
  unsigned long maxtry() const { return theMaxTry; }

  /**
   * Return the number of accepted points after which the grid should
   * be frozen
   */
  unsigned long freezeGrid() const { return theFreezeGrid; }

  /**
   * Set the number of points to presample this
   * splitting generator.
   */
  void presamplingPoints(unsigned long p) { thePresamplingPoints = p; }

  /**
   * Set the maximum number of trials
   * to generate a splitting.
   */
  void maxtry(unsigned long p) { theMaxTry = p; }

  /**
   * Set the number of accepted points after which the grid should
   * be frozen
   */
  void freezeGrid(unsigned long n) { theFreezeGrid = n; }

  /**
   * Evalute the splitting kernel.
   */
  double evaluate(const vector<double>&);

  /**
   * Return the index of the random number corresponding
   * to the evolution variable.
   */
  int evolutionVariable() const {
    return
      nDimBorn() +
      (showerApproximation()->showerInvertedTildeKinematics() ?
       showerApproximation()->showerInvertedTildeKinematics()->evolutionVariable() :
       dipole()->invertedTildeKinematics()->evolutionVariable());
  }

  /**
   * Return the cutoff on the evolution
   * random number corresponding to the pt cut.
   */
  double evolutionCutoff() const {
    return 
      showerApproximation()->showerInvertedTildeKinematics() ?
      showerApproximation()->showerInvertedTildeKinematics()->evolutionCutoff() :
      dipole()->invertedTildeKinematics()->evolutionCutoff();
  }

  /**
   * True, if sampler should apply compensation
   */
  void doCompensate(bool yes = true) { theDoCompensate = yes; }

public:

  /**@name Wrap to the exsample2 interface until this is finally cleaned up. */
  //@{

  inline const vector<bool>& variable_flags () {
    return sampleFlags();
  }

  inline size_t evolution_variable () const { return evolutionVariable(); }

  inline double evolution_cutoff () const { 
    return evolutionCutoff();
  }

  inline const vector<double>& parameter_point () {
    return parameterPoint();
  }

  inline void start_presampling () { 
    startPresampling();
  }

  inline void stop_presampling () { 
    stopPresampling();
  }

  inline size_t dimension () const { 
    return nDim();
  }

  inline unsigned long presampling_points () const { 
    return presamplingPoints();
  }

  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The dipole in charge of the emission
   */
  Ptr<SubtractionDipole>::ptr theDipole;

  /**
   * The shower approximation to consider
   */
  Ptr<ShowerApproximation>::ptr theShowerApproximation;

  /**
   * The XComb off which radiation will be generated
   */
  StdXCombPtr theBornXComb;

  /**
   * The XComb describing the process after radiation
   */
  StdXCombPtr theRealXComb;

  /**
   * The tilde xcomb objects associated to the real xcomb
   */
  vector<StdXCombPtr> theTildeXCombs;

  /**
   * True, if currently being presampled
   */
  bool thePresampling;

  /**
   * The number of points to presample this
   * splitting generator.
   */
  unsigned long thePresamplingPoints;

  /**
   * The maximum number of trials
   * to generate a splitting.
   */
  unsigned long theMaxTry;

  /**
   * Return the number of accepted points after which the grid should
   * be frozen
   */
  unsigned long theFreezeGrid;

  /**
   * The sampling flags
   */
  vector<bool> theFlags;

  /**
   * The support.
   */
  pair<vector<double>,vector<double> > theSupport;

  /**
   * The shower approximation generator.
   */
  Ptr<ShowerApproximationGenerator>::tptr theShowerApproximationGenerator;

  /**
   * The last parameter point
   */
  vector<double> theLastParameterPoint;

  /**
   * The last random numbers used for Born sampling
   */
  vector<double> theLastBornPoint;

  /**
   * Define the Sudakov sampler
   */
  typedef
  exsample::exponential_generator<ShowerApproximationKernel,UseRandom>
  ExponentialGenerator;

  /**
   * Define a pointer to the Sudakov sampler
   */
  typedef
  exsample::exponential_generator<ShowerApproximationKernel,UseRandom>*
  ExponentialGeneratorPtr;

  /**
   * The Sudakov sampler
   */
  ExponentialGeneratorPtr sampler;

  /**
   * True, if sampler should apply compensation
   */
  bool theDoCompensate;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerApproximationKernel & operator=(const ShowerApproximationKernel &);

};

}

#endif /* Herwig_ShowerApproximationKernel_H */
