// -*- C++ -*-
#ifndef HERWIG_TwoBodyAllOnCalculator_H
#define HERWIG_TwoBodyAllOnCalculator_H
// This is the declaration of the TwoBodyAllOnCalculator class.

#include "WidthCalculatorBase.h"
#include "TwoBodyAllOnCalculator.fh"
#include "GenericWidthGenerator.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>TwoBodyAllOnCalculator</code> class is a wrapped around the the 
 *  simple two body decay matrix elements in the <code>GenericWidthGenerator</code>
 *  class and is designed to allow these matrix elements to be integrated
 *  if the external particles can be off-shell.
 *
 * @see TwoBodyAllOnCalculator
 * 
 */
class TwoBodyAllOnCalculator: public WidthCalculatorBase {

public:

  /**
   * The GenericWidthGenerator class is a friend to allow easier access for the
   * integration of the two body partial widths.
   */
  friend class GenericWidthGenerator;

public:


  /** @name Standard constructors and destructors. */
  //@{
  /*
   * Constructor.
   * @param inwidth Pointer to the  GenericWidthGenerator class.
   * @param imode The mode in the GenericWidthGenerator class we are integrating
   * @param m1 The mass of the first particle.
   * @param m2 The mass of the second particle.
   */
  inline TwoBodyAllOnCalculator(tGenericWidthGeneratorPtr inwidth,int imode,
				Energy m1,Energy m2);

  /**
   * Default constructor
   */
  inline TwoBodyAllOnCalculator();

  /**
   * Copy constructor
   */
  inline TwoBodyAllOnCalculator(const TwoBodyAllOnCalculator &);

  /**
   * Destructor
   */
  virtual ~TwoBodyAllOnCalculator();
  //@}


public:

  /**
   * member to calculate the partial width.
   * @param scale The mass squared for the decaying particle.
   * @return The partial width.
   */
  Energy partialWidth(Energy2 scale) const;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @param mass The new value.
   * @return The mass required.
   */
  inline void resetMass(int imass,Energy mass);

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  inline Energy getMass(const int imass) const;

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  inline Energy otherMass(const int imass) const;

public:

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<TwoBodyAllOnCalculator> initTwoBodyAllOnCalculator;

  /**
   * Private and non-existent assignment operator.
   */
  TwoBodyAllOnCalculator & operator=(const TwoBodyAllOnCalculator &);

private:

  /**
   * the mode
   */
  int _mode;

  /**
   * Mass of the first particle.
   */
  Energy _mass1;

  /**
   * Mass of the second particle.
   */
  Energy _mass2;

  /**
   * the width generator
   */
  GenericWidthGeneratorPtr _widthgen;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of TwoBodyAllOnCalculator.
 */
template <>
 struct BaseClassTrait<Herwig::TwoBodyAllOnCalculator,1> {
  /** Typedef of the base class of TwoBodyAllShellCalculator. */
  typedef Herwig::WidthCalculatorBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TwoBodyAllOnCalculator>
  : public ClassTraitsBase<Herwig::TwoBodyAllOnCalculator> {
  /** Return the class name.*/
  static string className() { return "/Herwig++/TwoBodyAllOnCalculator"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return ""; }

};

}

#include "TwoBodyAllOnCalculator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoBodyAllOnCalculator.tcc"
#endif

#endif /* HERWIG_TwoBodyAllOnCalculator_H */
