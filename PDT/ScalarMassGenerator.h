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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ScalarMassGenerator();

  /**
   * Copy constructor
   */
  inline ScalarMassGenerator(const ScalarMassGenerator &);

  /**
   * Destructor
   */
  virtual ~ScalarMassGenerator();

public:

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentInput(PersistentIStream &, int);
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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
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
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return ""; }

};

}

#include "ScalarMassGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMassGenerator.tcc"
#endif

#endif /* HERWIG_ScalarMassGenerator_H */
