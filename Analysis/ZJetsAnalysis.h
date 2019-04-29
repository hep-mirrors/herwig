// -*- C++ -*-
#ifndef Herwig_ZJetsAnalysis_H
#define Herwig_ZJetsAnalysis_H
//
// This is the declaration of the ZJetsAnalysis class.
//

#include "Herwig/Analysis/JetsPlusAnalysis.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the ZJetsAnalysis class.
 *
 * @see \ref ZJetsAnalysisInterfaces "The interfaces"
 * defined for ZJetsAnalysis.
 */
class ZJetsAnalysis: public Herwig::JetsPlusAnalysis {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ZJetsAnalysis();

  /**
   * The destructor.
   */
  virtual ~ZJetsAnalysis();
  //@}

public:

  /**
   * Reconstruct the desired electroweak objects and fill the
   * respective momenta. Remove the reconstructed particles from the
   * list.
   */
  virtual void reconstructHardObjects(ParticleVector&);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ZJetsAnalysis & operator=(const ZJetsAnalysis &) = delete;

};

}

#endif /* Herwig_ZJetsAnalysis_H */
