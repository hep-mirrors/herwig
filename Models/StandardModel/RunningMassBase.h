// -*- C++ -*-
#ifndef HERWIG_RunningMassBase_H
#define HERWIG_RunningMassBase_H
//
// This is the declaration of the RunningMassBase class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Models
 * 
 *  Base class for running mass calculations.
 */
class RunningMassBase: public Interfaced {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline RunningMassBase();

  /**
   * Copy-constructor.
   */
  inline RunningMassBase(const RunningMassBase &);

  /**
   * Destructor.
   */
  virtual ~RunningMassBase();
  //@}
  
public:
  
  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param q2 The scale \f$q^2\f$.
   * @param part The ParticleData pointer
   */
  virtual Energy value(Energy2 q2,tcPDPtr part) const = 0;
 
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const = 0;

  /**
   * Return the \f$i\f$th element of the mass array.
   * @param i The element to return
   */
  inline Energy massElement(unsigned int i) const;

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
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<RunningMassBase> initRunningMassBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMassBase & operator=(const RunningMassBase &);
  
private:
  
  /**
   * Flavour thresholds and the masses, set at initialization.
   */
  vector<Energy> _theMass;

};

}

namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of RunningMassBase.
   */
  template <>
  struct BaseClassTrait<Herwig::RunningMassBase,1> {
    /** Typedef of the base class of RunningMassBase. */
    typedef Interfaced NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::RunningMassBase>
    : public ClassTraitsBase<Herwig::RunningMassBase> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::RunningMassBase"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwStandardModel.so"; }

  };
  
}

#include "RunningMassBase.icc"

#endif /* HERWIG_RunningMassBase_H */
