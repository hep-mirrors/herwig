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
  
  /**
   * Standard ctors and dtor.
   */
  inline RunningMass();
  inline RunningMass(const RunningMass &);
  virtual ~RunningMass();
  
public:
  
  /**
   * Return the running mass for a given scale and particle type.
   */
  virtual Energy value(Energy2 ,tcPDPtr) const;
  
  /**
   * Return the masses used.
   */
  virtual vector<Energy> mass() const;

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
   * Standard clone methods.
   */
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  
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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RunningMass> initRunningMass;
  
  /**
   * Private and non-existent assignment operator.
   */
  RunningMass & operator=(const RunningMass &);
  
  /**
   * Order in alphaS.
   */
  unsigned int _theQCDOrder;

  /**
   * The maximum number of flavours.
   */
  unsigned int _theMaxFlav;
  vector<double> _thePower,_theCoefficient;

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
    static string className() { return "/Herwig++/RunningMass"; }

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
