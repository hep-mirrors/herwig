// -*- C++ -*-
#ifndef HERWIG_RunningMass_H
#define HERWIG_RunningMass_H
//
// This is the declaration of the RunningMass class.

#include "RunningMassBase.h"
#include "StandardModel.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Models
 * 
 *  Implementation of the 1 or 2 loop QCD running mass.
 *
 *  @see RunningMassBase
 */
class RunningMass: public RunningMassBase {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline RunningMass();

  /**
   * Copy-constructor.
   */
  inline RunningMass(const RunningMass &);

  /**
   * Destructor.
   */
  virtual ~RunningMass();
  //@}  

public:
  
  /**
   * Return the running mass for a given scale \f$q^2\f$ and particle type.
   * @param q2 The scale \f$q^2\f$.
   * @param part The ParticleData pointer
   */
  virtual Energy value(Energy2 q2,tcPDPtr part) const;
  
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const;

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
  static ClassDescription<RunningMass> initRunningMass;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMass & operator=(const RunningMass &);
  
private:

  /**
   * Order in alphaS.
   */
  unsigned int _theQCDOrder;

  /**
   * The maximum number of flavours.
   */
  unsigned int _theMaxFlav;

  /**
   * The power for the running mass calculation.
   */
  vector<double> _thePower;

  /**
   * The coefficients for the running mass calculation.
   */
  vector<double> _theCoefficient;

  /**
   * Pointer to the StandardModel object.
   */
  SMPtr _theStandardModel;

};

}


namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of RunningMass.
   */
  template <>
  struct BaseClassTrait<Herwig::RunningMass,1> {
    /** Typedef of the base class of RunningMass. */
    typedef Herwig::RunningMassBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::RunningMass>
    : public ClassTraitsBase<Herwig::RunningMass> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::RunningMass"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwStandardModel.so"; }

  };
  
}

#include "RunningMass.icc"

#endif /* HERWIG_RunningMass_H */
