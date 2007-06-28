// -*- C++ -*-
#ifndef HERWIG_MSSM_H
#define HERWIG_MSSM_H
//
// This is the declaration of the MSSM class.
//

#include "SusyBase.h"
#include "MSSM.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MSSM class.
 *
 * @see \ref MSSMInterfaces "The interfaces"
 * defined for MSSM.
 */
class MSSM: public SusyBase {

public:

  /**
   * The default constructor.
   */
  inline MSSM();

public:

  /**
   * Value of Higgs mixing angle \f$\alpha\f$.
   */
  inline double higgsMixingAngle() const;

  /**
   * Value of up-type trilinear couplings
   */
  inline const complex<Energy> & topTrilinear() const;

  /**
   * Value of down-type trilinear couplings
   */
  inline const complex<Energy> & bottomTrilinear() const;

  /**
   * Value of lepton trilinear couplings
   */
  inline const complex<Energy> & tauTrilinear() const;

  /**
   * The stop mixing matrix
   */
  inline const MixingMatrixPtr & stopMix() const;

  /**
   * The sbottom chargino mixing matrix
   */
  inline const MixingMatrixPtr & sbottomMix() const;

  /**
   * The stau mixing matrix
   */
  inline const MixingMatrixPtr & stauMix() const;

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
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   * Adjust row of Mixing Matrix if a negative mass occurs in LHA file
   * @param id The PDG code of the particle with a negative mass
   */
  virtual void adjustMixingMatrix(long id);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MSSM> initMSSM;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MSSM & operator=(const MSSM &);

private:

  /**
   *  Third generation squark and slepton mixing matrices
   */
  //@{
  /**
   *  The \f$\tilde{t}\f$ mixing matrix
   */
  MixingMatrixPtr theStopMix;

  /**
   *  The \f$\tilde{b}\f$ mixing matrix
   */
  MixingMatrixPtr theSbotMix;

  /**
   *  The \f$\tilde{\tau}\f$ mixing matrix
   */
  MixingMatrixPtr theStauMix;
  //@}

  /**
   * Trilinear couplings stored as vector of complex numbers to make use
   * of routine already available to read complex matrices
   */
  //@{
  /**
   *  For the up type squarks
   */
  complex<Energy> theAtop;

  /**
   *  For the down type squarks
   */
  complex<Energy> theAbottom;

  /**
   *  For the charged sleptons
   */
  complex<Energy> theAtau;
  //@}

  /**
   * Value of higgs mixing angle.
   */
  double theAlpha;

  /**
   *  Higgs mixing matrix
   */
  MixingMatrixPtr theHiggsMix;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MSSM. */
template <>
struct BaseClassTrait<Herwig::MSSM,1> {
  /** Typedef of the first base class of MSSM. */
  typedef Herwig::SusyBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MSSM class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MSSM>
  : public ClassTraitsBase<Herwig::MSSM> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MSSM"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MSSM is implemented. It may also include several, space-separated,
   * libraries if the class MSSM depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "MSSM.icc"

#endif /* HERWIG_MSSM_H */
