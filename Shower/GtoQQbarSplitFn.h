// -*- C++ -*-
#ifndef HERWIG_GtoQQbarSplitFun_H
#define HERWIG_GtoQQbarSplitFun_H
//
// This is the declaration of the GtoQQbarSplitFun class.

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This abstract class provides the exact Leading Order splitting
 *  function for g -> q qbar . 
 *
 *  @see SplitFun1to2
 *  @see IS_GtoQQbarSplitFun
 *  @see FS_GtoQQbarSplitFun
 */
class GtoQQbarSplitFn: public SplittingFunction {

public:

  /**
   * Standard ctors and dtor.
   */
  inline GtoQQbarSplitFn() : SplittingFunction(ShowerIndex::QCD) {}
  inline GtoQQbarSplitFn(const GtoQQbarSplitFn &x) : SplittingFunction(x) {}
  virtual ~GtoQQbarSplitFn();

  /**
   * These virtual methods return the exact values of the 
   * Leading Order splitting function g -> q qbar
   * evaluated in terms of some combinations of:
   * z variable, phi azimuthal angle, and helicities of the three particles.
   */
  virtual double P(const double z, const Energy2 qtilde2, const IdList &);

  virtual double overestimateP(const double z, const IdList &);
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  /**
   * See long comment on this method on class SplitFun1to2.
   */
  virtual void colourConnection(const ShoColinePair &parent,
				ShoColinePair &first,
				ShoColinePair &second);

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

  static ClassDescription<GtoQQbarSplitFn> initGtoQQbarSplitFn; 

  /**
   * Private and non-existent assignment operator.
   */
  GtoQQbarSplitFn & operator=(const GtoQQbarSplitFn &);

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.
 */
template <>
struct BaseClassTrait<Herwig::GtoQQbarSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::GtoQQbarSplitFn>
  : public ClassTraitsBase<Herwig::GtoQQbarSplitFn>
{
  /** 
   * Return the class name.
   */
  static string className() { return "/Herwig++/GtoQQbarSplitFn"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwShower.so"; }

};

}

//#include "GtoQQbarSplitFun.icc"

#endif /* HERWIG_GtoQQbarSplitFn_H */
