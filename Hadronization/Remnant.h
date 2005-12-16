// -*- C++ -*-
#ifndef HERWIG_Remnant_H
#define HERWIG_Remnant_H
//
// This is the declaration of the Remnant class.
//

#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/PDF/PartonBinInstance.h>
#include <ThePEG/EventRecord/Step.h>
#include <ThePEG/EventRecord/Particle.h>
#include <Herwig++/Utilities/EnumParticles.h>
#include "Remnant.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The Remnant class is designed to store the hadronic remnant
 */
class Remnant: public Particle {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Remnant();

  /**
   *  Constructor from a parton bin instance
   * @param pb The parton bin instance
   * @param p The momentum
   */
  Remnant(PartonBinInstance & pb,const LorentzMomentum & p);

  /**
   *  Constructor with particle data object
   */
  Remnant(tcEventPDPtr);

  /**
   * The copy constructor.
   */
  inline Remnant(const Remnant &);

  /**
   * The destructor.
   */
  virtual ~Remnant();

  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void * operator new(size_t);
  
  /**
   * Particle uses the FixedSizeAllocator for (de)allocation.
   */
  inline void operator delete(void *, size_t);
  //@}

  /**
   *  Work out what the constituents
   *  @param extracted The PDG code of the extracted parton
   */
  void obtainConstituents(int extracted);

  /**
   *  Change the remnant if the extracted parton has change, after for
   *  example the initial-state shower
   * @param extracted The extracted parton
   * @param ptotal New momentum of the remnant
   */
  void regenerate(tPPtr extracted,Lorentz5Momentum ptotal);

  /**
   *  Make the constituents
   */
  void makeConstituents();

  /**
   *   create the remnant
   */
  void createRemnant(tStepPtr pstep);

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual PPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual PPtr fullclone() const;
  //@}

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static ClassDescription<Remnant> initRemnant;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Remnant & operator=(const Remnant &);

private:

  /**
   * Scale of the last perturbative branching
   */
  Energy2 _startscale;

  /**
   *  The incoming beam particle for which this is the remnant
   */
  cPDPtr _parent;

  /**
   *  The PDG code of the extracted particle
   */
  int _extracted;

  /**
   *  Whether this is the remnant for a particle or antiparticle
   */
  int _sign;

  /**
   *  The consituents of the remnant
   */
  ParticleVector _constituents;

  /**
   *  The id's of the valence partons
   */
  vector<int> _valence;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Remnant. */
template <>
struct BaseClassTrait<Herwig::Remnant,1> {
  /** Typedef of the first base class of Remnant. */
  typedef EventRecordBase  NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Remnant class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Remnant>
  : public ClassTraitsBase<Herwig::Remnant> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::Remnant"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the Remnant class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwMRST.so"; }
  /** Create a Event object. */
  static TPtr create() { return TPtr::Create(Herwig::Remnant(tcEventPDPtr())); }
};


/** @endcond */

}

#include "Remnant.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Remnant.tcc"
#endif

#endif /* HERWIG_Remnant_H */
