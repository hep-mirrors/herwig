// -*- C++ -*-
#ifndef HERWIG_QtoQGSplitFn_H
#define HERWIG_QtoQGSplitFn_H
//
// This is the declaration of the <!id>QtoQGSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting <BR>
// function for <I>Q-&GT;QG</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_QtoQGSplitFun.html">IS_QtoQGSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGSplitFun.html">FS_QtoQGSplitFun.h</a>.
// 

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

class QtoQGSplitFn: public SplittingFunction {

public:

  inline QtoQGSplitFn(const QtoQGSplitFn &x) : SplittingFunction(x) {}
  virtual ~QtoQGSplitFn();
  // Standard ctors and dtor.

  inline QtoQGSplitFn() : SplittingFunction(ShowerIndex::QCD) {}

  virtual double P(const double, const Energy2, const IdList &);

  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>Q-&GT;QG</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.

  virtual double overestimateP(const double z, const IdList &);
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  virtual void colourConnection(const ShoColinePair &parent,
				ShoColinePair &first, ShoColinePair &second);
  // See long comment on this method on class <!class>SplitFn1to2<!!class>.
  // Remember that the first branching product is considered the quark
  // and the second one the gluon.

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
  static ClassDescription<QtoQGSplitFn> initQtoQGSplitFn;  
  QtoQGSplitFn & operator=(const QtoQGSplitFn &);
  //  Private and non-existent assignment operator.

};

}

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ShowerHandler.
template <>
struct BaseClassTrait<Herwig::QtoQGSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::QtoQGSplitFn>
  : public ClassTraitsBase<Herwig::QtoQGSplitFn>
{
  static string className() { return "/Herwig++/QtoQGSplitFn"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}
//#include "QtoQGSplitFn.icc"

#endif /* HERWIG_QtoQGSplitFn_H */
