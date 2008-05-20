// -*- C++ -*-
#ifndef HERWIG_BtoCharmoniumDecayer_H
#define HERWIG_BtoCharmoniumDecayer_H
//
// This is the declaration of the BtoCharmoniumDecayer class.
//

#include "PartonicDecayerBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The BtoCharmoniumDecayer class is a simple model of the inclusive
 * decay of a b-hadron to a charmonium resonance and a hadronic system
 * which is hadronized using the cluster model.
 *
 * The decay is model using a spectator model approach where the b quark
 * decays to a charmonium resonance and a light quark, normally the strange quark.
 *
 *
 * @see \ref BtoCharmoniumDecayerInterfaces "The interfaces"
 * defined for BtoCharmoniumDecayer.
 */
class BtoCharmoniumDecayer: public PartonicDecayerBase {

public:

  /**
   * The default constructor.
   */
  BtoCharmoniumDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const PDVector & children) const;

  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const PDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BtoCharmoniumDecayer> initBtoCharmoniumDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoCharmoniumDecayer & operator=(const BtoCharmoniumDecayer &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BtoCharmoniumDecayer. */
template <>
struct BaseClassTrait<Herwig::BtoCharmoniumDecayer,1> {
  /** Typedef of the first base class of BtoCharmoniumDecayer. */
  typedef Herwig::PartonicDecayerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoCharmoniumDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoCharmoniumDecayer>
  : public ClassTraitsBase<Herwig::BtoCharmoniumDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BtoCharmoniumDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * BtoCharmoniumDecayer is implemented. It may also include several, space-separated,
   * libraries if the class BtoCharmoniumDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPartonicDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BtoCharmoniumDecayer_H */
