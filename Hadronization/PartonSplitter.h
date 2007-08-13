// -*- C++ -*-
#ifndef HERWIG_PartonSplitter_H
#define HERWIG_PartonSplitter_H

#include "CluHadConfig.h"
#include <ThePEG/Interface/Interfaced.h>
#include "PartonSplitter.fh"

namespace Herwig {


using namespace ThePEG;


/** \ingroup Hadronization
 *  \class PartonSplitter
 *  \brief This class splits the gluons from the end of the shower.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 * 
 *  This class does all of the nonperturbative parton splittings needed 
 *  immediately after the end of the showering (both initial and final),
 *  as very first step of the cluster hadronization.
 *
 *  \todo change so quark weights can be varied and quarks other
 *        than u and d can be produced
 */
class PartonSplitter: public Interfaced {

public:

  /**
   * This method does the nonperturbative splitting of:
   * time-like gluons. At the end of the shower the gluons should be
   * on a "physical" mass shell and should therefore be time-like.
   * @param tagged The tagged particles to be split
   * @return The particles which were not split and the products of splitting.
   */
  void split(PVector & tagged);
 
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<PartonSplitter> initPartonSplitter;

  /**
   * Private and non-existent assignment operator.
   */
  PartonSplitter & operator=(const PartonSplitter &);

  /**
   * Non-perturbatively split a time-like gluon,
   * if something goes wrong null pointers are returned.
   * @param gluon The gluon to be split
   * @param quark The quark produced in the splitting
   * @param anti  The antiquark produced in the splitting
   */
  void splitTimeLikeGluon(tcPPtr gluon, PPtr & quark, PPtr & anti);

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of PartonSplitter.
 */
template <>
struct BaseClassTrait<Herwig::PartonSplitter,1> {
  /** Typedef of the first base class of PartonSplitter. */
  typedef Interfaced NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::PartonSplitter>:
    public ClassTraitsBase<Herwig::PartonSplitter> {
  /** Return the class name.*/
  static string className() { return "Herwig::PartonSplitter"; }
};

/** @endcond */

}

#include "PartonSplitter.icc"

#endif /* HERWIG_PartonSplitter_H */
