// -*- C++ -*-
#ifndef HERWIG_GtoGGSplitFun_H
#define HERWIG_GtoGGSplitFun_H
//
// This is the declaration of the <!id>GtoGGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting <BR>
// function for <I>G-&GT;GG</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_GtoGGSplitFun.html">IS_GtoGGSplitFun.h</a>, <BR>
// <a href="http:FS_GtoGGSplitFun.html">FS_GtoGGSplitFun.h</a>.
//  

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

class GtoGGSplitFn: public SplittingFunction {

public:

  inline GtoGGSplitFn() : SplittingFunction(ShowerIndex::QCD) {}
  inline GtoGGSplitFn(const GtoGGSplitFn &x) : SplittingFunction(x) {}
  virtual ~GtoGGSplitFn();
  // Standard ctors and dtor.

  virtual double P(const double, const Energy2, const IdList &);
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>G-&GT;GG</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and
  // helicities of the three particles.

  virtual double overestimateP(const double z, const IdList &);
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );
  // See long comment on this method on class <!class>SplitFun1to2<!!class>.
  // Now interfaced routines
  static inline void Init() {}
  inline void persistentOutput(PersistentOStream &) const {}
  inline void persistentInput(PersistentIStream &, int) {} 
protected:
  inline virtual IBPtr clone() const { return new_ptr(*this); }
  inline virtual IBPtr fullclone() const { return clone(); }

  inline virtual void 
        doupdate() throw(UpdateException) { SplittingFunction::doupdate(); }
  inline virtual void 
        doinit() throw(InitException) { SplittingFunction::doinit(); }
  inline virtual void dofinish() { SplittingFunction::dofinish(); }
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException) { SplittingFunction::rebind(trans); }
  inline virtual IVector 
        getReferences() { return SplittingFunction::getReferences(); }
private:
  static ClassDescription<GtoGGSplitFn> initGtoGGSplitFn; 
  GtoGGSplitFn & operator=(const GtoGGSplitFn &);
  //  Private and non-existent assignment operator.

};

}

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ShowerHandler.
template <>
struct BaseClassTrait<Herwig::GtoGGSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::GtoGGSplitFn>
  : public ClassTraitsBase<Herwig::GtoGGSplitFn>
{
  static string className() { return "/Herwig++/GtoGGSplitFn"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}
#endif /* HERWIG_GtoGGSplitFun_H */
