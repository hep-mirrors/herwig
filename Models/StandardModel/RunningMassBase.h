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
  
  /**
   * Standard ctors and dtor.
   */
  inline RunningMassBase();
  inline RunningMassBase(const RunningMassBase &);
  virtual ~RunningMassBase();
  
public:
  
  /**
   * Return the running mass for a given scale and particle type.
   */
  virtual Energy value(Energy2 ,tcPDPtr) const = 0;
 
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const = 0;

  /**
   * Return the ith element of the mass array.
   */
  inline Energy massElement(unsigned int) const;

  /**
   * Standard functions for writing and reading from persistent streams.
   */
  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
protected:
  
  /**
   * Standard Interfaced virtual functions.
   */
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  
  /**
   * Change all pointers to Interfaced objects to corresponding clones.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  
  /**
   * Return pointers to all Interfaced objects refered to by this.
   */
  inline virtual IVector getReferences();
  
private:
  
  /**
   * Flavour thresholds and the masses, set at initialization.
   */
  vector<Energy> _theMass;

private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<RunningMassBase> initRunningMassBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMassBase & operator=(const RunningMassBase &);

};

}


namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of RunningMassBase.
   */
  template <>
  struct BaseClassTrait<Herwig::RunningMassBase,1> {
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
    static string className() { return "/Herwig++/RunningMassBase"; }

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
