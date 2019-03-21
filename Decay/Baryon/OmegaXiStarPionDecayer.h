// -*- C++ -*-
#ifndef HERWIG_OmegaXiStarPionDecayer_H
#define HERWIG_OmegaXiStarPionDecayer_H
// This is the declaration of the OmegaXiStarPionDecayer class.

#include "Baryon1MesonDecayerBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The OmegaXiStarPionDecayer class implements the results of
 *  hep-ph/0405162 for the weak decay of the \f$\Omega\f$ to the \f$\Xi^*\f$ and a pion.
 *
 * @see Baryon1MesonDecayerBase
 * 
 */
class OmegaXiStarPionDecayer: public Baryon1MesonDecayerBase {

public:

  /**
   * Default constructor.
   */
  OmegaXiStarPionDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /**
   *  Coupling Members.
   */
  //@{
  /**
   * Couplings for spin-\f$\frac32\f$ to spin-\f$\frac32\f$ and a scalar.
   * @param imode The mode
   * @param m0 The mass of the decaying particle.
   * @param m1 The mass of the outgoing baryon.
   * @param m2 The mass of the outgoing meson.
   * @param A1 The coupling \f$A_1\f$ described above.
   * @param A2 The coupling \f$A_2\f$ described above.
   * @param B1 The coupling \f$B_1\f$ described above.
   * @param B2 The coupling \f$B_2\f$ described above.
   */
  virtual void threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
						Complex& A1,Complex& A2,
						Complex& B1,Complex& B2) const;
  //@}

public:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<OmegaXiStarPionDecayer> initOmegaXiStarPionDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  OmegaXiStarPionDecayer & operator=(const OmegaXiStarPionDecayer &) = delete;

private:

  /**
   * The \f$A_{\rm Comm}\f$ amplitude from hep-ph/0405162
   */
  double Acomm_;

  /**
   * The \f$A_P\f$ amplitude from hep-ph/0405162
   */
  double AP_;

  /**
   * The \f$A_S\f$ amplitude from hep-ph/0405162
   */
  double AS_;

  /**
   * The \f$B_P\f$ amplitude from hep-ph/0405162
   */
  double BP_;

  /**
   * The \f$B_S\f$ amplitude from hep-ph/0405162
   */
  double BS_;
  
  /**
   * PDG code of the incoming baryon
   */
  int idin_;

  /**
   * PDG code of the outgoing baryon
   */
  int idout_;

  /**
   * maximum weight for the decay
   */
  double wgtmax_;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of OmegaXiStarPionDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::OmegaXiStarPionDecayer,1> {
    /** Typedef of the base class of OmegaXiStarPionDecayer. */
   typedef Herwig::Baryon1MesonDecayerBase NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
 struct ClassTraits<Herwig::OmegaXiStarPionDecayer>
  : public ClassTraitsBase<Herwig::OmegaXiStarPionDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::OmegaXiStarPionDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwBaryonDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_OmegaXiStarPionDecayer_H */
