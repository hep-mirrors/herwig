// -*- C++ -*-
#ifndef HERWIG_QtoGQSplitFn_H
#define HERWIG_QtoGQSplitFn_H
//
// This is the declaration of the QtoGQSplitFun class.

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This abstract class provides the exact Leading Order splitting
 *  function for <I>Q-&GT;QG</I>. 
 *
 *  @see SplitFun1to2
 *  @see IS_QtoGQSplitFun
 *  @see FS_QtoGQSplitFun
 */
class QtoGQSplitFn: public SplittingFunction {

public:

  /**
   * Standard ctors and dtor.
   */
  inline QtoGQSplitFn(const QtoGQSplitFn &x) : SplittingFunction(x) {}
  virtual ~QtoGQSplitFn();

  inline QtoGQSplitFn() : SplittingFunction(ShowerIndex::QCD) {}

  virtual double P(const double, const Energy2, const IdList &);

  /**
   * These virtual methods return the exact values of the 
   * Leading Order splitting function <I>Q-&GT;QG</I>
   * evaluated in terms of some combinations of:
   * z variable, phi azimuthal angle, and helicities of the three particles.
   */
  virtual double overestimateP(const double z, const IdList &);
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  /**
   * See long comment on this method on class  SplitFn1to2 .
   * Remember that the first branching product is considered the quark
   * and the second one the gluon.
   */
  virtual void colourConnection(const ShoColinePair &parent,
				ShoColinePair &first, ShoColinePair &second);

  /**
   * Now interfaced routines.
   */
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

  /**
   * Private and non-existent assignment operator.
   */
  static ClassDescription<QtoGQSplitFn> initQtoGQSplitFn;  
  QtoGQSplitFn & operator=(const QtoGQSplitFn &);

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.
 */
template <>
struct BaseClassTrait<Herwig::QtoGQSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::QtoGQSplitFn>
  : public ClassTraitsBase<Herwig::QtoGQSplitFn>
{
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/QtoGQSplitFn"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }

};

}
//#include "QtoGQSplitFn.icc"

#endif /* HERWIG_QtoGQSplitFn_H */
