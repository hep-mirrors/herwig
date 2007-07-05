// -*- C++ -*-
#ifndef HERWIG_NMSSM_H
#define HERWIG_NMSSM_H
//
// This is the declaration of the NMSSM class.
//

#include "MSSM.h"
#include "NMSSM.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the NMSSM class.
 *
 * @see \ref NMSSMInterfaces "The interfaces"
 * defined for NMSSM.
 */
class NMSSM: public MSSM {

public:

  /**
   * The default constructor.
   */
  inline NMSSM();

public:

  /**
   * Mixing matrix for the neutral CP-odd Higgs bosons
   */
  inline const MixingMatrixPtr & CPoddHiggsMix() const;

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  inline double lambda() const;

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  inline double kappa() const;
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

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<NMSSM> initNMSSM;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NMSSM & operator=(const NMSSM &);

private:

  /**
   *  Higgs mixing matrix
   */
  MixingMatrixPtr theHiggsAMix;

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  double _lambda;

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  double _kappa;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NMSSM. */
template <>
struct BaseClassTrait<Herwig::NMSSM,1> {
  /** Typedef of the first base class of NMSSM. */
  typedef Herwig::MSSM NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NMSSM class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NMSSM>
  : public ClassTraitsBase<Herwig::NMSSM> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::NMSSM"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NMSSM is implemented. It may also include several, space-separated,
   * libraries if the class NMSSM depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "NMSSM.icc"

#endif /* HERWIG_NMSSM_H */
