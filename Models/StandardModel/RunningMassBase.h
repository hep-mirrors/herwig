// -*- C++ -*-
#ifndef HERWIG_RunningMassBase_H
#define HERWIG_RunningMassBase_H
//
// This is the declaration of the <!id>RunningMassBase<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
//  Base class for running mass calculations
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:.html">.h</a>,
// <a href="http:.html">.h</a>.
// 

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

class RunningMassBase: public Interfaced {
  
public:
  
  inline RunningMassBase();
  inline RunningMassBase(const RunningMassBase &);
  virtual ~RunningMassBase();
  // Standard ctors and dtor.
  
public:
  
  virtual Energy value(Energy2 ,tcPDPtr) const = 0;
  // Return the running mass for a given scale and particle type.

  virtual vector<Energy> mass() const = 0;
  // Return the masses used.

  inline Energy massElement(unsigned int) const;
  // return the ith element of the mass array


  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
protected:
  
  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void doinitrun();
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.
  
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.
  
  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.
  
private:
  
  vector<Energy> _theMass;
  // flavour thresholds and the masses, set at initialization.

private:
  
  static AbstractClassDescription<RunningMassBase> initRunningMassBase;
  // Describe an abstract base class with persistent data.
  
  RunningMassBase & operator=(const RunningMassBase &);
  // Private and non-existent assignment operator.

};

}

// CLASSDOC OFF

namespace ThePEG {

  // The following template specialization informs ThePEG about the
  // base class of RunningMassBase.
  template <>
  struct BaseClassTrait<Herwig::RunningMassBase,1> {
    typedef Interfaced NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::RunningMassBase>
    : public ClassTraitsBase<Herwig::RunningMassBase> {
    static string className() { return "/Herwig++/RunningMassBase"; }
    // Return the class name.
    static string library() { return "libHwStandardModel.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "RunningMassBase.icc"

#endif /* HERWIG_RunningMassBase_H */
