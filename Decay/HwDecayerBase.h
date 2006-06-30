// -*- C++ -*-
#ifndef HERWIG_HwDecayerBase_H
#define HERWIG_HwDecayerBase_H
//
// This is the declaration of the HwDecayerBase class.
//

#include "ThePEG/PDT/Decayer.h"
#include "HwDecayerBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The HwDecayerBase class is the base class for Decayers in Herwig++. It inherits
 * from the Decayer class of ThePEG and implements additional functionality for the 
 * output of the results to the particle database and initialization of the datbase.
 *
 * It also provide the option of specifying a class based on the DecayRadiationGenerator
 * which should be used to generate QED radiation in the decay
 *
 */
class HwDecayerBase: public Decayer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline HwDecayerBase();
  //@}

public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:

  /**
   *  Functions for the Herwig decayer
   */
  //@{

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

  /**
   *  Access to the initialize variable
   */
  inline bool initialize() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<HwDecayerBase> initHwDecayerBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwDecayerBase & operator=(const HwDecayerBase &);

private:

  /**
   * perform initialisation
   */
  bool _initialize;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HwDecayerBase. */
template <>
struct BaseClassTrait<Herwig::HwDecayerBase,1> {
  /** Typedef of the first base class of HwDecayerBase. */
  typedef Decayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HwDecayerBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HwDecayerBase>
  : public ClassTraitsBase<Herwig::HwDecayerBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::HwDecayerBase"; }
  /** Return the name of the shared library be loaded to get
   *  access to the HwDecayerBase class and every other class it uses
   *  (except the base class). */
  static string library() { return ""; }
};

/** @endcond */

}

#include "HwDecayerBase.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HwDecayerBase.tcc"
#endif

#endif /* HERWIG_HwDecayerBase_H */
