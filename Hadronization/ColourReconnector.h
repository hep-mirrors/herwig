// -*- C++ -*-
#ifndef HERWIG_ColourReconnector_H
#define HERWIG_ColourReconnector_H

#include <ThePEG/Handlers/HandlerBase.h>
#include "CluHadConfig.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Hadronization
 *  \class ColourReconnector
 *  \brief Class for changing colour reconnections of partons.
 *  \author Alberto Ribon
 * 
 *  This class does the nonperturbative colour rearrangement, after the 
 *  nonperturbative gluon splitting and the "normal" cluster formation. 
 *  It uses the list of particles in the event record, and the collections of
 *  "usual" clusters which is passed to the main method. If the colour 
 *  reconnection is actually accepted, then the previous collections of "usual"
 *  clusters is first deleted and then the new one is created.
 *
 *  Note: by default this class does nothing. It can be inherited and overridden
 *  in future hadronization models.
 */
//class ThePEG::PartialCollisionHandler; // forward declaration


class ColourReconnector: public ThePEG::HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline ColourReconnector();

  /**
   * Copy-constructor.
   */
  inline ColourReconnector(const ColourReconnector &);

  /**
   * Destructor.
   */
  virtual ~ColourReconnector();
  //@}

  /**
   * Does the colour rearrangment.
   *
   * Does the colour rearrangement, starting from the list of particles
   * in the event record, and the collection of "usual" clusters passed
   * in input. If the actual rearrangement is accepted, the new collection 
   * of clusters is overriden to the intial one.
   */
  void rearrange(EventHandler & ch, const StepPtr & pstep,
                 ClusterVector & clusters) throw(Veto, Stop, Exception);
    
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
  static ClassDescription<ColourReconnector> initColourReconnector;

  /**
   * Private and non-existent assignment operator.
   */
  ColourReconnector & operator=(const ColourReconnector &);

  /**
   * The numer of colours in the reconstruction.
   */
  int    _ClReco;

  /**
   * The probability of a reconstruction.
   */
  double _PReco;

};


}

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of ColourReconnector.
 */
struct BaseClassTrait<Herwig::ColourReconnector,1> {
  /** Typedef of the base class of ColourReconnector. */
  typedef ThePEG::HandlerBase NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::ColourReconnector>
  : public ClassTraitsBase<Herwig::ColourReconnector> {
  /** Return the class name.*/
  static string className() { return "Herwig++::ColourReconnector"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwHadronization.so"; }

};

}

#endif // DOXYGEN

#include "ColourReconnector.icc"

#endif /* HERWIG_ColourReconnector_H */
