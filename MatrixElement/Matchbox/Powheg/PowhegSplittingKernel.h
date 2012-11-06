// -*- C++ -*-
//
// PowhegSplittingKernel.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegSplittingKernel_H
#define HERWIG_PowhegSplittingKernel_H
//
// This is the declaration of the PowhegSplittingKernel class.
//

#include "Herwig++/MatrixElement/Matchbox/Powheg/ME2byDipoles.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"

namespace Herwig {

using namespace ThePEG;

class PowhegSplittingGenerator;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief PowhegSplittingKernel implements the splitting
 * kernel entering the POWHEG Sudakov form factor.
 *
 */
class PowhegSplittingKernel: public ME2byDipoles {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PowhegSplittingKernel();

  /**
   * The destructor.
   */
  virtual ~PowhegSplittingKernel();
  //@}

public:

  /**
   * Evaluate the ratio.
   */
  virtual double evaluate() const;

  /**
   * Return true, if 'Born screening' should be done
   */
  bool bornScreening() const { return theBornScreening; }

  /**
   * Switch on 'Born screening'
   */
  void doBornScreening() { theBornScreening = true; }

  /**
   * Switch off 'Born screening'
   */
  void noBornScreening() { theBornScreening = false; }

  /**
   * Set a pt cut
   */
  void ptCut(Energy pt) { projectionDipole()->ptCut(pt); }

  /**
   * Set a screening scale
   */
  void screeningScale(Energy s) { theScreeningScale = s; }

  /**
   * Set a reference to the splitting generator
   * making use of this kernel.
   */
  void splittingGenerator(Ptr<PowhegSplittingGenerator>::tptr gen);

public:

  /**
   * Set the XComb object steering the real emission.
   */
  virtual void setXComb(tStdXCombPtr real);

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches();

  /**
   * Return the number of random numbers
   * needed to sample this kernel.
   */
  int nDim() const;

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
   * Evalute the splitting kernel.
   */
  double evaluate(const vector<double>&);

  /**
   * Return true, if currently being presampled
   */
  bool presampling() const { return thePresampling; }

  /**
   * Return the index of the random number corresponding
   * to the evolution variable.
   */
  int evolutionVariable() const;

  /**
   * Return the cutoff on the evolution
   * random number corresponding to the pt cut.
   */
  double evolutionCutoff() const {
    return projectionDipole()->invertedTildeKinematics()->evolutionCutoff();
  }

public:

  /**
   * Return the subprocess from which splittings
   * have been generated.
   */
  tSubProPtr bornSubProcess() const {
    return projectionDipole()->lastHeadXCombPtr()->subProcess();
  }

  /**
   * Return a new subprocess, if selected for splitting
   */
  tSubProPtr construct(Energy pt);

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
   * True, if 'Born screening' should be done
   */
  bool theBornScreening;

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
   * An optional screening scale; see
   * DipoleSplittingKernel
   */
  Energy theScreeningScale;

  /**
   * The range of born random numbers
   */
  pair<int,int> theBornRandom;

  /**
   * The range of radiation random numbers
   */
  pair<int,int> theRadiationRandom;

  /**
   * The last parameter point
   */
  vector<double> theLastParameterPoint;

  /**
   * Random numbers used to presample the Born
   */
  vector<double> thePresamplingPoint;

  /**
   * True, if currently being presampled
   */
  bool thePresampling;

  /**
   * Backup of the Born XComb when being presampled.
   */
  StdXCombPtr theXCombBackup;

  /**
   * Map Born XCombs to private ones used for presampling.
   */
  map<StdXCombPtr,StdXCombPtr> thePresamplingXCombs;

  /**
   * A reference to the splitting generator
   * making use of this kernel.
   */
  Ptr<PowhegSplittingGenerator>::tptr theGenerator;

  /**
   * The sampling flags
   */
  vector<bool> theFlags;

  /**
   * The support.
   */
  pair<vector<double>,vector<double> > theSupport;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegSplittingKernel & operator=(const PowhegSplittingKernel &);

};

}

#endif /* HERWIG_PowhegSplittingKernel_H */
