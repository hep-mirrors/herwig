// -*- C++ -*-
#ifndef HERWIG_QtoQGammaSplitFun_H
#define HERWIG_QtoQGammaSplitFun_H
//
// This is the declaration of the <!id>QtoQGammaSplitFun<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This abstract class provides the exact Leading Order splitting
// function for <I>Q-&GT;QGamma</I>. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:SplitFun1to2.html">SplitFun1to2.h</a>, <BR>
// <a href="http:IS_QtoQGammaSplitFun.html">IS_QtoQGammaSplitFun.h</a>, <BR>
// <a href="http:FS_QtoQGammaSplitFun.html">FS_QtoQGammaSplitFun.h</a>.
//  

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

class QtoQGammaSplitFn: public SplittingFunction {

public:

  inline QtoQGammaSplitFn() : SplittingFunction(ShowerIndex::QED) {}
  inline QtoQGammaSplitFn(const QtoQGammaSplitFn &x) : SplittingFunction(x) {}
  virtual ~QtoQGammaSplitFn();
  // Standard ctors and dtor.

  virtual double P(const double z, const Energy2 qtilde2, const IdList &);
  // These virtual methods return the exact values of the 
  // Leading Order splitting function <I>Q-&GT;QGamma</I>
  // evaluated in terms of some combinations of:
  // <!id>z<!!id> variable, <!id>phi<!!id> azimuthal angle, and 
  // helicities of the three particles.

  virtual double overestimateP(const double z, const IdList &); 
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  virtual void colourConnection(const ShoColinePair &parent,
				ShoColinePair &first,
				ShoColinePair &second);
  // See long comment on this method on class <!class>SplitFun1to2<!!class>.
  // Remember that the first branching product is considered the quark
  // and the second one the photon.
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
  static ClassDescription<QtoQGammaSplitFn> initQtoQGammaSplitFn; 
  QtoQGammaSplitFn & operator=(const QtoQGammaSplitFn &);
  //  Private and non-existent assignment operator.

};

}


namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ShowerHandler.
template <>
struct BaseClassTrait<Herwig::QtoQGammaSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::QtoQGammaSplitFn>
  : public ClassTraitsBase<Herwig::QtoQGammaSplitFn>
{
  static string className() { return "/Herwig++/QtoQGammaSplitFn"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}
//#include "QtoQGammaSplitFun.icc"

#endif /* HERWIG_QtoQGammaSplitFun_H */
