// -*- C++ -*-
#ifndef HERWIG_ExampleDecayer_H
#define HERWIG_ExampleDecayer_H
//
// This is the declaration of the ExampleDecayer class.
//
#include <ThePEG/PDT/Decayer.h>

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  Dummy class which provides just an example of how to inherits
   *  from the abstract class  Decayer.  and then implement the
   *  physics of the decay in the two methods:
   *            accept
   *            decay
   *
   * 
   */
class ExampleDecayer: public ThePEG::Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ExampleDecayer();

  /**
   * Copy constructor
   */
  inline ExampleDecayer(const ExampleDecayer &);

  /**
   * Destructor
   */
  virtual ~ExampleDecayer();
  //@}

public:

  /**
   * return true if this decayer can perfom the decay specified by the
   * given decay mode.
   */
  virtual bool accept(const DecayMode &) const;

  /**
   * for a given decay mode and a given particle instance, perform the
   * decay and return the decay products.
   */
  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
 
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
  static ClassDescription<ExampleDecayer> initExampleDecayer;

  /**
   *  Private and non-existent assignment operator.
   */
  ExampleDecayer & operator=(const ExampleDecayer &);

};


}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ExampleDecayer.
 */
template <>
struct BaseClassTrait<Herwig::ExampleDecayer,1> {
  /** Typedef of the base class of ExampleDecayer. */
  typedef ThePEG::Decayer NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ExampleDecayer>: public ClassTraitsBase<Herwig::ExampleDecayer> {
  /** Return the class name. */
  static string className() { return "Herwig++::ExampleDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwDecay.so"; }

};

}

#include "ExampleDecayer.icc"

#endif /* HERWIG_ExampleDecayer_H */
