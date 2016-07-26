// -*- C++ -*-
#ifndef Herwig_PowhegHandler_H
#define Herwig_PowhegHandler_H
//
// This is the declaration of the PowhegHandler class.
//

#include "MatchingHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the PowhegHandler class.
 *
 * @see \ref PowhegHandlerInterfaces "The interfaces"
 * defined for PowhegHandler.
 */
class PowhegHandler: public MatchingHandler {

public:

  /**
   * The default constructor.
   */
  PowhegHandler()  : MatchingHandler(false), pTDefinition_(0), maxpT_(GeV)
  {}

  /**
   * Perform CKKW reweighting
   */
  virtual double reweightCKKW(int minMult, int maxMult);

  /**
   *  Generate hard emissions for CKKW etc
   */
  virtual HardTreePtr generateCKKW(ShowerTreePtr tree) const;

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
   *  Calculate the Sudakov weight
   */
  virtual double sudakovWeight( CKKWTreePtr ) {
    return 1.;
  }

  virtual PotentialTree chooseHardTree(double totalWeight,
				       double nonOrderedWeight);

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
  PowhegHandler & operator=(const PowhegHandler &);

private:

  /**
   *  Control over the selection of the maximum \f$p_T\f$ for emission
   */
  unsigned int pTDefinition_;

  /**
   *  Maximum pT for the shower for non-emission events
   */
  Energy maxpT_;

};

}

#endif /* Herwig_PowhegHandler_H */
