// -*- C++ -*-
#ifndef HERWIG_RunningMass_H
#define HERWIG_RunningMass_H
//
// This is the declaration of the <!id>RunningMass<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// Implementation of the 1 or 2 loop QCD running mass 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="RunningMassBase.html">RunningMaseBase.h</a>.
// 

#include "RunningMassBase.h"
#include "StandardModel.h"

namespace Herwig {
using namespace ThePEG;

class RunningMass: public RunningMassBase {
  
public:
  
  inline RunningMass();
  inline RunningMass(const RunningMass &);
  virtual ~RunningMass();
  // Standard ctors and dtor.
  
public:
  
  virtual Energy value(Energy2 ,tcPDPtr) const;
  // Return the running mass for a given scale and particle type.
  
  virtual vector<Energy> mass() const;
  // Return the masses used.

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.
  
  static void Init();
  // Standard Init function used to initialize the interfaces.
  
protected:
  
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.
  
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
  
  static ClassDescription<RunningMass> initRunningMass;
  // Describe a concrete class with persistent data.
  
  RunningMass & operator=(const RunningMass &);
  // Private and non-existent assignment operator.
  
  unsigned int _theQCDOrder;
  // order in alphaS
  unsigned int _theMaxFlav;
  vector<double> _thePower,_theCoefficient;
  // the maximum number of flavours
  SMPtr _theStandardModel;
  // pointer to the StandardModel object
};

}

// CLASSDOC OFF

namespace ThePEG {
  
  // The following template specialization informs ThePEG about the
  // base class of RunningMass.
  template <>
  struct BaseClassTrait<Herwig::RunningMass,1> {
    typedef Herwig::RunningMassBase NthBase;
  };
  
  // The following template specialization informs ThePEG about the
  // name of this class and the shared object where it is defined.
  template <>
  struct ClassTraits<Herwig::RunningMass>
    : public ClassTraitsBase<Herwig::RunningMass> {
    static string className() { return "/Herwig++/RunningMass"; }
    // Return the class name.
    static string library() { return "libHwStandardModel.so"; }
    // Return the name of the shared library to be loaded to get
    // access to this class and every other class it uses
    // (except the base class).
  };
  
}

#include "RunningMass.icc"

#endif /* HERWIG_RunningMass_H */
