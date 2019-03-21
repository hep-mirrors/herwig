// -*- C++ -*-
//
// DipoleMCCheck.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleMCCheck_H
#define HERWIG_DipoleMCCheck_H
//
// This is the declaration of the DipoleMCCheck class.
//

#include "ThePEG/Handlers/HandlerBase.h"

#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 * 
 * \brief DipoleMCCheck is used to perform checks for
 * the dipole shower.
 *
 * @see \ref DipoleMCCheckInterfaces "The interfaces"
 * defined for DipoleMCCheck.
 */
class DipoleMCCheck: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleMCCheck();

  /**
   * The destructor.
   */
  virtual ~DipoleMCCheck();
  //@}

public:

  /**
   * Book a point.
   */
  void book(double xe,double xs,
	    Energy dScale,
	    Energy hardPt,
	    Energy pt, double z,
	    double weight);

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



protected:

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The number of bins in the starting scale
   * divided by the dipole scale, the upper bound
   * is 1/2.
   */
  unsigned int theHardPtBins;

  /**
   * The number of bins in the emitter fraction.
   * The lenght of the zero bin is taken to be
   * 10^(-7).
   */
  unsigned int theEmitterXBins;

  /**
   * The number of bins in the spectator fraction.
   * The lenght of the zero bin is taken to be
   * 10^(-7).
   */
  unsigned int theSpectatorXBins;

  /**
   * The number of bins in pt dicided by the
   * dipole scale; the upper bound is 1/2
   */
  unsigned int thePtBins;

  /**
   * The number of bins in z
   */
  unsigned int theZBins;

  /**
   * The recursive map structure: xe, xs, hard pt / GeV
   * to histograms for pt and z;
   * output is done such that there's one file for each
   * xe,xs bin, including all histograms for the 
   * hard pt bins.
   */
  map<double,
      map<double,
	  map<double,
	      pair<Ptr<Histogram>::ptr,Ptr<Histogram>::ptr>
	      >
	  >
      > histoMap;

  /**
   * Helper to make logarithmic bins.
   */
  vector<double> makeLogBins(double xlow, double xup, unsigned int n) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleMCCheck> initDipoleMCCheck;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleMCCheck & operator=(const DipoleMCCheck &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleMCCheck. */
template <>
struct BaseClassTrait<Herwig::DipoleMCCheck,1> {
  /** Typedef of the first base class of DipoleMCCheck. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleMCCheck class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleMCCheck>
  : public ClassTraitsBase<Herwig::DipoleMCCheck> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleMCCheck"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleMCCheck is implemented. It may also include several, space-separated,
   * libraries if the class DipoleMCCheck depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleMCCheck_H */
