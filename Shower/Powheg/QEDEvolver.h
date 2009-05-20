// -*- C++ -*-
#ifndef HERWIG_QEDEvolver_H
#define HERWIG_QEDEvolver_H
//
// This is the declaration of the QEDEvolver class.
//

#include "Herwig++/Shower/Base/Evolver.h"

namespace Herwig {

using namespace ThePEG;

  /**\ingroup Shower
   * Exception class
   * used to communicate failure of QED shower
   */
  struct QEDVeto {};

/**
 * Here is the documentation of the QEDEvolver class.
 *
 * @see \ref QEDEvolverInterfaces "The interfaces"
 * defined for QEDEvolver.
 */
class QEDEvolver: public Evolver {

public:

  /**
   * The default constructor.
   */
  QEDEvolver();

  /**
   *  Members to perform the shower
   */
  //@{
  /**
   * Perform the shower of the hard process
   */
  virtual void showerHardProcess(ShowerTreePtr);

  /**
   * Perform the shower of a decay
   */
  virtual void showerDecay(ShowerTreePtr);
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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
   *   Construct a time-like line
   */
  void constructTimeLikeLine(tHardBranchingPtr branch,tShowerParticlePtr particle);

  void constructSpaceLikeLine(tShowerParticlePtr particle,
			      HardBranchingPtr & first, HardBranchingPtr & last,
			      SudakovPtr sud,PPtr beam);

  void constructHardTree(vector<ShowerProgenitorPtr> & particlesToShower,
			 ShowerInteraction::Type inter);
  
  bool constructDecayTree(vector<ShowerProgenitorPtr> & particlesToShower,
			  ShowerInteraction::Type inter);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<QEDEvolver> initQEDEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QEDEvolver & operator=(const QEDEvolver &);

private:

  /**
   *  Interactions allowed in the shower
   */
  vector<ShowerInteraction::Type> interactions_;

  /**
   *  Order of the interactions
   */
  bool QCDFirst_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QEDEvolver. */
template <>
struct BaseClassTrait<Herwig::QEDEvolver,1> {
  /** Typedef of the first base class of QEDEvolver. */
  typedef Herwig::Evolver NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QEDEvolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QEDEvolver>
  : public ClassTraitsBase<Herwig::QEDEvolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QEDEvolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QEDEvolver is implemented. It may also include several, space-separated,
   * libraries if the class QEDEvolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QEDEvolver_H */
