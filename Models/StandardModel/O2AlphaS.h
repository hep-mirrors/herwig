// -*- C++ -*-
#ifndef HERWIG_O2AlphaS_H
#define HERWIG_O2AlphaS_H
//
// This is the declaration of the O2AlphaS class.
//

#include "ThePEG/StandardModel/AlphaSBase.h"
#include "O2AlphaS.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The O2AlphaS class is the implementation of the two-loop
 * \f$\alpha_S\f$ in the same way as in FORTRAN HERWIG.
 *
 *  The input value of \f$\Lambda_{\rm QCD}\f$ is in the \f$\bar{MS}\f$
 *  scheme and can either be converted to a Monte Carlo scheme, as in the 
 *  FORTRAN program, or left in the \f$\bar{MS}\f$ scheme to evaluate
 *  the running coupling 
 *
 * @see \ref O2AlphaSInterfaces "The interfaces"
 * defined for O2AlphaS.
 */
class O2AlphaS: public AlphaSBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline O2AlphaS();

  /**
   * The copy constructor.
   */
  inline O2AlphaS(const O2AlphaS &);

  /**
   * The destructor.
   */
  virtual ~O2AlphaS();
  //@}

public:

  /** @name Virtual functions to override those in the base class */
  //@{
  /**
   * The \f$\alpha_S\f$. Return the QCD coupling for a given \a scale
   * using the given standard model object \a sm.
   */
  virtual double value(Energy2 scale, const StandardModelBase & sm) const;

  /**
   * Return the flavour thresholds used. The returned vector contains
   * (in position <code>i</code>) the scales when the active number of
   * flavours changes from <code>i</code> to <code>i+1</code>.
   */
  virtual vector<Energy2> flavourThresholds() const;

  /**
   * Return the \f$\Lambda_{QCD}\f$ used for different numbers of
   * active flavours.
   */
  virtual vector<Energy> LambdaQCDs() const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);
  //@}


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<O2AlphaS> initO2AlphaS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  O2AlphaS & operator=(const O2AlphaS &);

private:

  /**
   *  The value of \f$\Lambda_{\rm QCD}\f$ 5-flavours
   *  in the \f$\bar{\rm MS}\f$ scheme
   */
  Energy _lambdaQCD;

  /**
   *  The values of the leading-order \f$\beta\f$-function coefficients
   */
  double _bcoeff[6];

  /**
   *  The values of the next-to-leading-order \f$\beta\f$-function coefficients
   */
  double _ccoeff[6];

  /**
   *  The values of \f$\Lambda_{\rm QCD}\f$ for the diffferent number of flavours
   */
  Energy _lambdas[6];

  /**
   *  The flavour thresholds
   */
  Energy _threshold[6];
  
  /**
   *  The constants for matching
   */
  double _match[6];
  /**
   *  Option for the coupling
   */
  unsigned int _copt;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of O2AlphaS. */
template <>
struct BaseClassTrait<Herwig::O2AlphaS,1> {
  /** Typedef of the first base class of O2AlphaS. */
  typedef AlphaSBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the O2AlphaS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::O2AlphaS>
  : public ClassTraitsBase<Herwig::O2AlphaS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::O2AlphaS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * O2AlphaS is implemented. It may also include several, space-separated,
   * libraries if the class O2AlphaS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwStandardModel.so"; }
};

/** @endcond */

}

#include "O2AlphaS.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "O2AlphaS.tcc"
#endif

#endif /* HERWIG_O2AlphaS_H */
