// -*- C++ -*-
#ifndef HERWIG_UEBase_H
#define HERWIG_UEBase_H
//
// This is the declaration of the UEBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardXComb.fh"
#include "UEBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Abstract base class used to minimize the dependence between
 * MPIHandler and all Shower classes.
 *
 * \author Manuel B\"ahr
 *
 * @see \ref UEBaseInterfaces "The interfaces"
 * defined for UEBase.
 */
class UEBase: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  UEBase(){}

  /**
   * The destructor.
   */
  virtual ~UEBase(){}
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /** @name Virtual functions used for the generation of additional
      interactions . */
  //@{

  /**
   * Some initialization code eventually.
   */
  virtual void initialize() {}

  /**
   * Return true or false depending on the generator setup. 
   */
  virtual bool beamOK() const = 0;

  /**
   * Return the value of the pt cutoff.
   */
  virtual Energy Ptmin() const = 0;

  /**
   * Return the slope of the soft pt spectrum. Only necessary when the
   * soft part is modelled.
   */
  virtual InvEnergy2 beta() const {return 0/GeV2;}

  /**
   * Some finalize code eventually.
   */
  virtual void finalize() {}

  /**
   * Return the number of different hard processes. Use 0 as default to
   * not require implementation.
   */
  virtual unsigned int additionalHardProcs() const {return 0;}

  /**
   * return the hard multiplicity of process i. Can't be constant in my
   * case because drawing from the probability distribution also
   * specifies the soft multiplicity that has to be stored....
   */
  virtual unsigned int multiplicity(unsigned int i=0) = 0;

  /**
   * Generate a additional interaction for ProcessHandler sel. Method
   * can't be const because it saves the state of the underlying XComb
   * object on it's way.
   */
  virtual tStdXCombPtr generate(unsigned int sel=0) = 0;

  /**
   * Return the type of algorithm. 
   */
  virtual int Algorithm() const = 0;

  /**
   * Return the value of the hard Process pt cutoff for vetoing.
   */
  virtual Energy PtForVeto() const = 0;

  /**
   * Return the fraction of colour disrupted subprocesses. Use default 0
   * so that it is not required to implement.
   */
  virtual double colourDisrupt() const {return 0.0;}

  /**
   * Return the soft multiplicity. Use 0 as default to not require
   * implementation.
   */
  virtual unsigned int softMultiplicity() const {return 0;} 

  //@}
  

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<UEBase> initUEBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEBase & operator=(const UEBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEBase. */
template <>
struct BaseClassTrait<Herwig::UEBase,1> {
  /** Typedef of the first base class of UEBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEBase>
  : public ClassTraitsBase<Herwig::UEBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEBase is implemented. It may also include several, space-separated,
   * libraries if the class UEBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "UEBase.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEBase_H */
