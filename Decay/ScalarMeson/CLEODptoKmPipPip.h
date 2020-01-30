// -*- C++ -*-
#ifndef Herwig_CLEODptoKmPipPip_H
#define Herwig_CLEODptoKmPipPip_H
//
// This is the declaration of the CLEODptoKmPipPip class.
//

#include "WeakDalitzDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The CLEODptoKmPipPip class implements the CLEO model from arXiv:0802.4214 for \f$D^+\toK^-\pi^+\pi^+\f$
 *
 * @see \ref CLEODptoKmPipPipInterfaces "The interfaces"
 * defined for CLEODptoKmPipPip.
 */
class CLEODptoKmPipPip: public WeakDalitzDecay {

public:

  /**
   * The default constructor.
   */
  CLEODptoKmPipPip();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  Calculate the amplitude
   */
  virtual Complex amplitude(int ichan) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CLEODptoKmPipPip & operator=(const CLEODptoKmPipPip &);

};

}

#endif /* Herwig_CLEODptoKmPipPip_H */
