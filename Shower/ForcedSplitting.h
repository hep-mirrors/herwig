// -*- C++ -*-
#ifndef HERWIG_ForcedSplitting_H
#define HERWIG_ForcedSplitting_H
//
// This is the declaration of the ForcedSplitting class.

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ThePEG/Repository/CurrentGenerator.h"
//#include "BackwardEvolver.h"

namespace Herwig {

class BackwardEvolver;

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This class is used to take the remaining parton from the backward evolution
 *  and tie it to the hadron with the appropriate remnant such that the
 *  last parton in the shower and the remnant constitute the complete set of
 *  flavours in the hadron 
 *
 *  @see BackwardEvolver
 */
class ForcedSplitting: public Interfaced {

public:

  /**
   * Standard ctors and dtor.
   */
  inline ForcedSplitting();
  inline ForcedSplitting(const ForcedSplitting &);
  virtual ~ForcedSplitting() {}
  int split(tShowerParticlePtr &part, ShowerParticleVector &p, tEHPtr ch);

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);

  static void Init();

private:

  static ClassDescription<ForcedSplitting> initForcedSplitting;
  ForcedSplitting & operator=(const ForcedSplitting &);


  Ptr<BackwardEvolver>::pointer _backEvolve;
  ShowerParticleVector *_particles;

  /**
   * This routine is used to generate the ShowerKinematics object in the 
   * forced splitting.
   */
  ShoKinPtr forcedSplitting(const ShowerParticle &, Energy, Energy);
  

  tEGPtr generator() const;
  Ptr<SplittingGenerator>::pointer splittingGenerator();
  
  // Sets the pointer so the proper list of particles are updated
  void setParticleList(ShowerParticleVector &);

  // This takes the particle and find a splitting for np -> p + child and 
  // creates the correct kinematics and connects for such a split. This
  // Splitting has an upper bound on qtilde given by the energy argument
  int forceSplit(tShowerParticlePtr &, long p, long np, long child, Energy &);

  // This creates a remnant for connected with the parton passed in. The
  // maxIdx is the number of partons in the incoming hadron and q[3] is
  // a vector which contains all of the flavours of the hadrons partons.
  // The last argument is the parton which is given by part.
  void makeRemnant(tShowerParticlePtr &part, int maxIdx, long q[3], int,
		   tEHPtr ch);

protected:
  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  inline virtual void rebind(const TranslationMap &) throw(RebindException);
  inline virtual IVector getReferences();
};

}


namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ShowerHandler.   
 */
template <>
struct BaseClassTrait<Herwig::ForcedSplitting,1> {
  typedef ThePEG::Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ForcedSplitting>
  : public ClassTraitsBase<Herwig::ForcedSplitting>
{

  /**
   * Return the class name.
   */
  static string className() { return "/Herwig++/ForcedSplitting"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwShower.so"; }
};

}

#include "ForcedSplitting.icc"

#endif /* HERWIG_ForcedSplitting_H */

