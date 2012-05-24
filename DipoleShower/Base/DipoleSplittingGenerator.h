// -*- C++ -*-
//
// DipoleSplittingGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleSplittingGenerator_H
#define HERWIG_DipoleSplittingGenerator_H
//
// This is the declaration of the DipoleSplittingGenerator class.
//

#include "ThePEG/Handlers/HandlerBase.h"

#include "Herwig++/DipoleShower/Kernels/DipoleSplittingKernel.h"
#include "DipoleSplittingReweight.h"
#include "Herwig++/DipoleShower/Utility/DipoleMCCheck.h"
#include "Herwig++/Exsample2/exsample/exponential_generator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleSplittingGenerator is used by the dipole shower
 * to sample splittings from a given dipole splitting kernel.
 *
 * @see \ref DipoleSplittingGeneratorInterfaces "The interfaces"
 * defined for DipoleSplittingGenerator.
 */
class DipoleSplittingGenerator: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleSplittingGenerator();

  /**
   * The destructor.
   */
  virtual ~DipoleSplittingGenerator();
  //@}

public:

  /**
   * Return the dipole splitting kernel.
   */
  Ptr<DipoleSplittingKernel>::tptr splittingKernel() const;

  /**
   * Return the dipole splitting reweight.
   */
  Ptr<DipoleSplittingReweight>::tptr splittingReweight() const;

  /**
   * Return the dipole splitting kinematics.
   */
  Ptr<DipoleSplittingKinematics>::tptr splittingKinematics() const;

  /**
   * Set the dipole splitting kernel.
   */
  void splittingKernel(Ptr<DipoleSplittingKernel>::tptr sp);

  /**
   * Set the dipole splitting reweight.
   */
  void splittingReweight(Ptr<DipoleSplittingReweight>::tptr sp);

  /**
   * Make a wrapper around another generator.
   */
  void wrap(Ptr<DipoleSplittingGenerator>::ptr other);

  /**
   * Return true, if this is actually a wrapper around
   * another splitting generator.
   */
  bool wrapping() const { return theOtherGenerator; }

public:

  /**
   * Prepare to fill the given splitting.
   */
  void prepare(const DipoleSplittingInfo&);

  /**
   * Fix parameters from the fiven DipoleSplittingInfo
   * and generate the next splitting. Return the
   * pt selected for the next splitting.
   */
  Energy generate(const DipoleSplittingInfo&);

  /**
   * Fix parameters from the fiven DipoleSplittingInfo
   * and generate the next splitting. Return the
   * pt selected for the next splitting when called
   * from a wrapping generator.
   */
  Energy generateWrapped(DipoleSplittingInfo&);

  /**
   * Complete the given splitting.
   */
  void completeSplitting(DipoleSplittingInfo&) const;

  /**
   * Return the last generated splitting
   */
  const DipoleSplittingInfo& lastSplitting() const { return generatedSplitting; }

public:

  /**
   * Print debug information on the splitting
   * handled.
   */
  void debugGenerator(ostream&) const;

  /**
   * Print debug information on the last
   * generated event.
   */
  void debugLastEvent(ostream&) const;

protected:

  /**
   * Update parameters given a splitting.
   */
  void fixParameters(const DipoleSplittingInfo&);

  /**
   * With the parameters previuosly supplied
   * through fixParameters generate the next
   * splitting.
   */
  void doGenerate();

public:

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
  const vector<double>& parameterPoint() const { return parameters; }

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
  unsigned long presamplingPoints() const { return splittingKernel()->presamplingPoints(); }

  /**
   * Return the maximum number of trials
   * to generate a splitting.
   */
  unsigned long maxtry() const { return splittingKernel()->maxtry(); }

  /**
   * Return true, if this splitting generator
   * is able to deliver an overestimate to the sampled
   * kernel.
   */
  bool haveOverestimate() const;

  /**
   * Return an overestimate to the sampled kernel.
   */
  bool overestimate(const vector<double>&);

  /**
   * Invert the integral over the overestimate to equal
   * the given value.
   */
  double invertOverestimateIntegral(double) const;

  /**
   * Evalute the splitting kernel.
   */
  double evaluate(const vector<double>&);

public:

  /**@name Wrap to the exsample2 interface until this is finally cleaned up. */
  //@{

  inline const vector<bool>& variable_flags () {
    return sampleFlags();
  }

  inline size_t evolution_variable () const { return 0; }

  inline double evolution_cutoff () { return support().first[0]; }

  inline const vector<double>& parameter_point () const {
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
   * Pointer to another generator to wrap around.
   */
  Ptr<DipoleSplittingGenerator>::ptr theOtherGenerator;

  /**
   * The dipole splitting kernel to sample
   * splitting from.
   */
  Ptr<DipoleSplittingKernel>::ptr theSplittingKernel;

  /**
   * The dipole splitting reweight.
   */
  Ptr<DipoleSplittingReweight>::ptr theSplittingReweight;

  /**
   * Pointer to the exponential generator
   */
  exsample::exponential_generator<DipoleSplittingGenerator,UseRandom>*
  theExponentialGenerator;

  /**
   * The dipole splitting to be completed.
   */
  DipoleSplittingInfo generatedSplitting;

  /**
   * A backup of the dipole splitting to be
   * completed, if this generator is presampled.
   */
  DipoleSplittingInfo presampledSplitting;

  /**
   * True, if prepared to sample splittings
   * of a given kind.
   */
  bool prepared;

  /**
   * Wether or not the kernel is currently
   * being presampled.
   */
  bool presampling;

  /**
   * The parameter point.
   */
  vector<double> parameters;

  /**
   * The sampling flags
   */
  vector<bool> theFlags;

  /**
   * The support.
   */
  pair<vector<double>,vector<double> > theSupport;

  /**
   * Pointer to a check histogram object
   */
  Ptr<DipoleMCCheck>::ptr theMCCheck;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleSplittingGenerator> initDipoleSplittingGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleSplittingGenerator & operator=(const DipoleSplittingGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleSplittingGenerator. */
template <>
struct BaseClassTrait<Herwig::DipoleSplittingGenerator,1> {
  /** Typedef of the first base class of DipoleSplittingGenerator. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleSplittingGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleSplittingGenerator>
  : public ClassTraitsBase<Herwig::DipoleSplittingGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleSplittingGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleSplittingGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DipoleSplittingGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleSplittingGenerator_H */
