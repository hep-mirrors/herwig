// -*- C++ -*-
#ifndef HERWIG_Hw64Decayer_H
#define HERWIG_Hw64Decayer_H
//
// This is the declaration of the <!id>Hw64Decayer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>Hw64Decayer<!!id> is a class that defines all the general routines 
// used in HERWIG++ to imitate the HERWIG 6.4 decays. The goal is to have an exact
// copy of HERWIG 6.4 decay routines. This will allow for easy 'callibration'
// of the new C++ code with the old Fortran code.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="Decayer.html">Decayer.h</a>.
// 

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/PDT/Decayer.h>
#include <ThePEG/Handlers/HandlerBase.h>
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/Strategy.fh>
#include "DecayConfig.h"
#include <fstream>

namespace Herwig {

using namespace ThePEG;

class Hw64Decayer: public Decayer {

public:

  inline Hw64Decayer();
  inline Hw64Decayer(const Hw64Decayer &);
  virtual ~Hw64Decayer();
  // Standard ctors and dtor

public:

  virtual bool accept(const DecayMode &) const;
  // return true if this decayer can perfom the decay specified by the
  // given decay mode.

  virtual ParticleVector decay(const DecayMode &, const Particle &) const;
  // for a given decay mode and a given particle instance, perform the
  // decay and return the decay products.

  static void Init();
  // Standard Init function used to initialize the interface.

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard Persisten stream methods

protected:
   inline virtual IBPtr clone() const;
   inline virtual IBPtr fullclone() const;
   // Standard clone methods

   inline virtual void rebind(const TranslationMap &trans) 
                              throw(RebindException);
   // Change all pointers of Interfaced objects to corresponding clones

   inline virtual IVector getReferences();
   // Return pointers to all interfaced objects referred to by this class.

private:

  int MECode;

  static void oneBodyDecay(Lorentz5Momentum, Lorentz5Momentum &);
  // Perform a one body decay, used for K
  // two body decay is handled in static class Kinematics

  static double VAWt(double*); //, double, double, double) const;
  // Weighting of phase space for V-A matrix elements

  //double PhaseSpaceWt() const;
  // Flat phase space weight (1.0)

  //double randomAzimuthal(double, double &, double &) const;
  void setParticleMomentum(ParticleVector &, ParticleMSet, 
			   vector<Lorentz5Momentum>) const;
  // Sets the momentum in the vector to be the momentum of the particle vector
  // where is the particles are of the type given in the ParticleMSet

  static ClassDescription<Hw64Decayer> initHw64Decayer;

  const Hw64Decayer & operator=(const Hw64Decayer &);
  //  Private and non-existent assignment operator.
};

}

namespace ThePEG {

template <>
struct BaseClassTrait<Herwig::Hw64Decayer,1> {
  typedef HandlerBase NthBase;
};

template <>
struct ClassTraits<Herwig::Hw64Decayer>: public ClassTraitsBase<Herwig::Hw64Decayer> {
  static string className() { return "/Herwig++/Hw64Decayer"; }
  static string library() { return "libHwDecay.so"; }
};

}

#include "Hw64Decayer.icc"

#endif /* HERWIG_Hw64Decayer_H */
