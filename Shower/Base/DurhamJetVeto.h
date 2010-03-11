// -*- C++ -*-
#ifndef HERWIG_DurhamJetVeto_H
#define HERWIG_DurhamJetVeto_H
//
// This is the declaration of the DurhamJetVeto class.
//

#include "ShowerVeto.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DurhamJetVeto class.
 *
 * @see \ref DurhamJetVetoInterfaces "The interfaces"
 * defined for DurhamJetVeto.
 */
class DurhamJetVeto: public ShowerVeto {

public:

  /**
   * The default constructor.
   */
  DurhamJetVeto() : ShowerVeto(ShowerVeto::Emission),
		    _ptVetoDefinition(1), _reversePtVeto(false), 
		    _y_cut(1.1), _approxCuts( false )  {}

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			     const Branching&);

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			      const Branching&);

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DurhamJetVeto> initDurhamJetVeto;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DurhamJetVeto & operator=(const DurhamJetVeto &);

private:

  /**
   * The pt definition being used for in the pt veto
   */
  unsigned int _ptVetoDefinition;

  /**
   * Reverse pt veto to produce branchings above a cut only
   */
  bool _reversePtVeto;
  
  /**
   *  \f$y_{\rm cut}\f$
   */
  double _y_cut;

  /**
   * switch to use approximate jet cuts
   */ 
  bool _approxCuts;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DurhamJetVeto. */
template <>
struct BaseClassTrait<Herwig::DurhamJetVeto,1> {
  /** Typedef of the first base class of DurhamJetVeto. */
  typedef Herwig::ShowerVeto NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DurhamJetVeto class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DurhamJetVeto>
  : public ClassTraitsBase<Herwig::DurhamJetVeto> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DurhamJetVeto"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DurhamJetVeto is implemented. It may also include several, space-separated,
   * libraries if the class DurhamJetVeto depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DurhamJetVeto_H */
