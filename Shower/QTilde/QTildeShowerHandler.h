// -*- C++ -*-
#ifndef Herwig_QTildeShowerHandler_H
#define Herwig_QTildeShowerHandler_H
//
// This is the declaration of the QTildeShowerHandler class.
//

#include "QTildeShowerHandler.fh"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Shower/QTilde/Base/Evolver.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The QTildeShowerHandler class.
 *
 * @see \ref QTildeShowerHandlerInterfaces "The interfaces"
 * defined for QTildeShowerHandler.
 */
class QTildeShowerHandler: public ShowerHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QTildeShowerHandler();

  /**
   * The destructor.
   */
  virtual ~QTildeShowerHandler();
  //@}

public:

  /**
   *  Access to the Evolver
   */
  tEvolverPtr evolver() const {return evolver_;}

  /**
   * At the end of the Showering, transform ShowerParticle objects
   * into ThePEG particles and fill the event record with them.
   * Notice that the parent/child relationships and the 
   * transformation from ShowerColourLine objects into ThePEG
   * ColourLine ones must be properly handled.
   */
  void fillEventRecord();
  
  /**
   * Return the relevant hard scale to be used in the profile scales
   */
  virtual Energy hardScale() const;

  /**
   *  Generate hard emissions for CKKW etc
   */
  virtual HardTreePtr generateCKKW(ShowerTreePtr tree) const;

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

  /**
   * The main method which manages the showering of a subprocess.
   */
  virtual tPPair cascade(tSubProPtr sub, XCPtr xcomb);

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeShowerHandler & operator=(const QTildeShowerHandler &);

private:

  /**
   *  Whether or not to split into hard and decay trees
   */
  bool splitHardProcess_;

  /**
   *  The ShowerTree for the hard process
   */
  ShowerTreePtr hard_;

  /**
   *  Pointer to the evolver
   */
  EvolverPtr evolver_;

  /**
   *  The ShowerTree for the decays
   */
  ShowerDecayMap decay_;

  /**
   *  The ShowerTrees for which the initial shower 
   */
  vector<ShowerTreePtr> done_;

};

}

#endif /* Herwig_QTildeShowerHandler_H */
