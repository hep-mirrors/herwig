// -*- C++ -*-
#ifndef HERWIG_GtoGGSplitFun_H
#define HERWIG_GtoGGSplitFun_H
//
// This is the declaration of the GtoGGSplitFun class.

#include "SplittingFunction.h"


namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * 
 *  This abstract class provides the exact Leading Order splitting
 *  function for g -> g g . 
 *
 *  @see SplitFun1to2
 *  @see IS_GtoGGSplitFun
 *  @see FS_GtoGGSplitFun
 */
class GtoGGSplitFn: public SplittingFunction {

public:

  /**
   * Standard ctors and dtor.
   */
  inline GtoGGSplitFn() : SplittingFunction(ShowerIndex::QCD) {}
  inline GtoGGSplitFn(const GtoGGSplitFn &x) : SplittingFunction(x) {}
  virtual ~GtoGGSplitFn();

  /**
   * These virtual methods return the exact values of the 
   * Leading Order splitting function g -> g g
   * evaluated in terms of some combinations of:
   * z variable, phi azimuthal angle, and helicities of the 
   * three particles.
   */
  virtual double P(const double, const Energy2, const IdList &);

  virtual double overestimateP(const double z, const IdList &);
  virtual double integOverP(const double z); 
  virtual double invIntegOverP(const double r); 

  /**
   * See long comment on this method on class SplitFun1to2.
   * Now interfaced routines.
   */
  virtual void colourConnection( const ShoColinePair & parentShoColinePair,
				 ShoColinePair & firstProductShoColinePair,
				 ShoColinePair & secondProductShoColinePair );

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
  static ClassDescription<GtoGGSplitFn> initGtoGGSplitFn; 
  GtoGGSplitFn & operator=(const GtoGGSplitFn &);

};

}

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.
 */
template <>
struct BaseClassTrait<Herwig::GtoGGSplitFn,1> {
  typedef Herwig::SplittingFunction NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::GtoGGSplitFn>
  : public ClassTraitsBase<Herwig::GtoGGSplitFn>
{
  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/GtoGGSplitFn"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }

};

}
#endif /* HERWIG_GtoGGSplitFun_H */
