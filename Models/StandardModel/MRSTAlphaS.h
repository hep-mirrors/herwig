// -*- C++ -*-
#ifndef HERWIG_MRSTAlphaS_H
#define HERWIG_MRSTAlphaS_H
//
// This is the declaration of the MRSTAlphaS class.
//

#include "ThePEG/StandardModel/AlphaSBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MRSTAlphaS class.
 *
 * @see \ref MRSTAlphaSInterfaces "The interfaces"
 * defined for MRSTAlphaS.
 */
class MRSTAlphaS: public AlphaSBase {

public:

  /**
   * The default constructor.
   */
  MRSTAlphaS();

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MRSTAlphaS> initMRSTAlphaS;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MRSTAlphaS & operator=(const MRSTAlphaS &);

private:

  /**
   * \f$\Lambda_{\rm QCD}\f$
   */
  Energy _lambda;

  /**
   * Charm Mass 
   */
  Energy _mcharm;

  /**
   * Bottom Mass
   */
  Energy _mbottom;

  /**
   *  Number of flavours
   */
  int _flavour;

  /**
   *  Order 
   */
  int _iord;

  /**
   *  \f$4m_c^2\f$
   */
  Energy2 _qsct;

  /**
   *  \f$4m_b^2\f$
   */
  Energy2 _qsbt;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MRSTAlphaS. */
template <>
struct BaseClassTrait<Herwig::MRSTAlphaS,1> {
  /** Typedef of the first base class of MRSTAlphaS. */
  typedef AlphaSBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MRSTAlphaS class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MRSTAlphaS>
  : public ClassTraitsBase<Herwig::MRSTAlphaS> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MRSTAlphaS"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MRSTAlphaS is implemented. It may also include several, space-separated,
   * libraries if the class MRSTAlphaS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMRSTAlphaS.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MRSTAlphaS_H */
