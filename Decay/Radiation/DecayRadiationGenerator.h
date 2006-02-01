// -*- C++ -*-
#ifndef HERWIG_DecayRadiationGenerator_H
#define HERWIG_DecayRadiationGenerator_H
//
// This is the declaration of the DecayRadiationGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "DecayRadiationGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The DecayRadiationGenerator class is the base class for classes generating
 * QED radiation in particle decays in Herwig++. Classes implementing specific
 * algorithms must inherit from this class and implement the virtual generatePhotons
 * member.
 */
class DecayRadiationGenerator: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DecayRadiationGenerator();

  /**
   * The copy constructor.
   */
  inline DecayRadiationGenerator(const DecayRadiationGenerator &);

  /**
   * The destructor.
   */
  virtual ~DecayRadiationGenerator();
  //@}

public:

  /**
   *  Member to generate the photons in the decay. This must be implemented
   *  in classes inheriting from this one to produce the radiation.
   * @param p The decaying particle
   * @param children The decay products
   * @return The decay products with additional radiation
   */
  virtual ParticleVector generatePhotons(const Particle & p,ParticleVector children)=0;

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
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
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
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DecayRadiationGenerator> initDecayRadiationGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DecayRadiationGenerator & operator=(const DecayRadiationGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DecayRadiationGenerator. */
template <>
struct BaseClassTrait<Herwig::DecayRadiationGenerator,1> {
  /** Typedef of the first base class of DecayRadiationGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DecayRadiationGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DecayRadiationGenerator>
  : public ClassTraitsBase<Herwig::DecayRadiationGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DecayRadiationGenerator"; }
  /** Return the name of the shared library be loaded to get
   *  access to the DecayRadiationGenerator class and every other class it uses
   *  (except the base class). */
  static string library() { return "DecayRadiationGenerator.so"; }
};

/** @endcond */

}

#include "DecayRadiationGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayRadiationGenerator.tcc"
#endif

#endif /* HERWIG_DecayRadiationGenerator_H */
