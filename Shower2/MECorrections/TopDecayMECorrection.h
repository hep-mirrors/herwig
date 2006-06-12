// -*- C++ -*-
#ifndef HERWIG_TopDecayMECorrection_H
#define HERWIG_TopDecayMECorrection_H
//
// This is the declaration of the TopDecayMECorrection class.
//

#include "MECorrectionBase.h"
#include "TopDecayMECorrection.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The TopDecayMECorrection class implements the matrix element correction
 * for top decay.
 *
 * @see \ref TopDecayMECorrectionInterfaces "The interfaces"
 * defined for TopDecayMECorrection.
 */
class TopDecayMECorrection: public MECorrectionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline TopDecayMECorrection();

  /**
   * The copy constructor.
   */
  inline TopDecayMECorrection(const TopDecayMECorrection &);

  /**
   * The destructor.
   */
  virtual ~TopDecayMECorrection();
  //@}

public:

  /**
   *  Members to override those in the base class and implemented 
   *  the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   */
  virtual bool canHandle(ShowerTreePtr);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br);
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /**
   *  Full matrix element with a factor of \f$\frac{\alpha_SC_F}{x_g^2\pi}\f$ removed.
   * @param xw The momentum fraction of the W boson
   * @param xg The momentum fraction of the gluon.
   */
  double me(double xw, double xg);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TopDecayMECorrection> initTopDecayMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TopDecayMECorrection & operator=(const TopDecayMECorrection &);

private:

  /**
   *  The mass of the W boson
   */
  Energy _mw;

  /**
   *  The mass of the bottom quark
   */
  Energy _mb;

  /**
   *  The top mass
   */
  Energy _mt;

  /**
   *  The mass ratio for the W
   */
  double _a;

  /**
   *  The mass ratio for the bottom
   */
  double _c;

  /**
   *  The enhancement factor for initial-state radiation
   */
  double _initialenhance;

  /**
   *  The enhancement factor for final-state radiation
   */
  double _finalenhance;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TopDecayMECorrection. */
template <>
struct BaseClassTrait<Herwig::TopDecayMECorrection,1> {
  /** Typedef of the first base class of TopDecayMECorrection. */
  typedef Herwig::MECorrectionBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TopDecayMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TopDecayMECorrection>
  : public ClassTraitsBase<Herwig::TopDecayMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::TopDecayMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TopDecayMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class TopDecayMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "TopDecayMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDecayMECorrection.tcc"
#endif

#endif /* HERWIG_TopDecayMECorrection_H */
