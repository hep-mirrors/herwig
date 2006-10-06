// -*- C++ -*-
#ifndef HERWIG_ScalarMassGenerator_H
#define HERWIG_ScalarMassGenerator_H
// This is the declaration of the ScalarMassGenerator class.

#include "GenericMassGenerator.h"
#include "ScalarMassGenerator.fh"
#include "ThePEG/Config/Complex.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup PDT
 *
 *  The <code>ScalarMassGenerator</code> class is designed for the generation
 *  of the masses of the \f$a_0\f$ and \f$f_0\f$ mesons which have \f$K\bar{K}\f$
 *  modes close
 *  to the on-shell mass of the particle. It includes finite-width effects.
 *
 * @see MassGenerator
 * @see GenericMassGenerator
 * 
 */
class ScalarMassGenerator: public GenericMassGenerator {

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

public:

  /**
   * Weight for the factor for an off-shell mass
   * @param mass The off-shell mass
   * @return The weight.
   */
  inline virtual double weight(Energy mass) const;

  /**
   * output for the database
   */
  virtual void dataBaseOutput(ofstream &);

protected:

  /**
   * The self-energy for the weight
   * @param mass The off-shell mass.
   * @return The self energy.
   */
  inline complex<Energy2> selfEnergy(Energy mass) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}


private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ScalarMassGenerator> initScalarMassGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarMassGenerator & operator=(const ScalarMassGenerator &);

 private:

  /**
   * couplings for the decay channels
   */
  vector<Energy> _coupling;

  /**
   * The first outgoing mass for the channels
   */
  vector<Energy> _mass1;

  /**
   * The second outgoing mass for the channels
   */
  vector<Energy> _mass2;


  /**
   * calculated values to speed things up
   */
  vector<Energy>  _mplus ,_mminus;

  /**
   * calculated values to speed things up
   */
  vector<Energy2> _m2plus,_m2minus;

  /**
   * calculated values to speed things up
   */
  vector<Energy2> _term;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ScalarMassGenerator.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarMassGenerator,1> {
  /** Typedef of the base class of ScalarMassGenerator. */
  typedef Herwig::GenericMassGenerator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ScalarMassGenerator>
  : public ClassTraitsBase<Herwig::ScalarMassGenerator> {
  /** Return the class name. */
  static string className() { return "/Herwig++/ScalarMassGenerator"; }
};

}

#include "ScalarMassGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMassGenerator.tcc"
#endif

#endif /* HERWIG_ScalarMassGenerator_H */
