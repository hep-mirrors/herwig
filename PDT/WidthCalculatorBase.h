// -*- C++ -*-
#ifndef HERWIG_WidthCalculatorBase_H
#define HERWIG_WidthCalculatorBase_H
//
// This is the declaration of the WidthCalculatorBase class.
//
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Interface/Interfaced.h"
#include "WidthCalculatorBase.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 * *

 *  The <code>WidthCalculatorBase</code> class is a base class to be used
 *  by classes which calculate partial widths for the running width.
 *
 * @see DecayIntegrator
 * @see GenericWidthGenerator
 * 
 */
class WidthCalculatorBase: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline WidthCalculatorBase();
  /**
   * Copy constructor
   */
  inline WidthCalculatorBase(const WidthCalculatorBase &);

  /**
   * Destructor
   */
  virtual ~WidthCalculatorBase();
  //@}

public:

  /**
   * Calculate the partial width. This must be implemented in classes inheriting from
   * this one.
   * @param scale The mass squared of the decaying particle.
   * @return The partial width.
   */
  virtual Energy partialWidth(Energy2 scale) const =0;

  /**
   * Reset the mass of a particle (used to integrate over the mass.) This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass to be reset.
   * @param mass The new mass.
   */
  virtual void resetMass(int imass,Energy mass) =0;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  virtual Energy getMass(const int imass) const= 0;

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  virtual Energy otherMass(const int imass) const=0;

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

private:

  /**
   * Describe an abstract base class without persistent data.
   */
  static AbstractNoPIOClassDescription<WidthCalculatorBase> initWidthCalculatorBase;

  /**
   * Private and non-existent assignment operator.
   */
  WidthCalculatorBase & operator=(const WidthCalculatorBase &);

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of WidthCalculatorBase.
 */
template <>
struct BaseClassTrait<Herwig::WidthCalculatorBase,1> {
  /** Typedef of the base class of WidthCalculatorBase. */
  typedef Base NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::WidthCalculatorBase>
  : public ClassTraitsBase<Herwig::WidthCalculatorBase> {
  /** Return the class name. */
  static string className() { return "Herwig++::WidthCalculatorBase"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return ""; }

};

}

#include "WidthCalculatorBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "WidthCalculatorBase.tcc"
#endif

#endif /* HERWIG_WidthCalculatorBase_H */
