// -*- C++ -*-
#ifndef HERWIG_MECorrectionBase_H
#define HERWIG_MECorrectionBase_H
//
// This is the declaration of the MECorrectionBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower2/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower2/ShowerConfig.h"
#include "Herwig++/Shower2/ShowerProgenitor.h"
#include "Herwig++/Shower2/ShowerTree.h"
#include "MECorrectionBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MECorrectionBase class.
 *
 * @see \ref MECorrectionBaseInterfaces "The interfaces"
 * defined for MECorrectionBase.
 */
class MECorrectionBase: public Interfaced {

/**
 *  The Evolver is a friend to set some variables at initialisation
 */
friend class Evolver;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MECorrectionBase();

  /**
   * The copy constructor.
   */
  inline MECorrectionBase(const MECorrectionBase &);

  /**
   * The destructor.
   */
  virtual ~MECorrectionBase();
  //@}

public:

  /**
   *  Virtual members which should be implemented in classes inheriting
   * from this one in order to apply the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   */
  virtual bool canHandle(ShowerTreePtr)=0;

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr)=0;

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br)=0;
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
   *  Set/Get methods for the pointer to ShowerVariables
   */
  //@{
  /**
   *  Get the ShowerVariables
   */
  inline void showerVariables(ShowerVarsPtr);
  
  /**
   *  Set the ShowerVariables
   */
  inline ShowerVarsPtr showerVariables() const;
  //@}  

  /**
   *  Access to the coupling
   */
  inline ShowerAlphaPtr coupling() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractClassDescription<MECorrectionBase> initMECorrectionBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MECorrectionBase & operator=(const MECorrectionBase &);

private:

  /**
   *  Pointer to the ShowerVariables object
   */
  ShowerVarsPtr _showerVariables;

  /**
   *  Pointer to the coupling
   */
  ShowerAlphaPtr _alpha;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MECorrectionBase. */
template <>
struct BaseClassTrait<Herwig::MECorrectionBase,1> {
  /** Typedef of the first base class of MECorrectionBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MECorrectionBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MECorrectionBase>
  : public ClassTraitsBase<Herwig::MECorrectionBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MECorrectionBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MECorrectionBase is implemented. It may also include several, space-separated,
   * libraries if the class MECorrectionBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNewShower.so"; }
};

/** @endcond */

}

#include "MECorrectionBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MECorrectionBase.tcc"
#endif

#endif /* HERWIG_MECorrectionBase_H */
