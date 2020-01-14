// -*- C++ -*-
#ifndef Herwig_KMatrix_H
#define Herwig_KMatrix_H
//
// This is the declaration of the KMatrix class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Decay/IsoSpin.h"
#include "KMatrix.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The KMatrix class is a base class for the implementation of
 * K-matrix parameterizations in Herwig
 *
 * @see \ref KMatrixInterfaces "The interfaces"
 * defined for KMatrix.
 */
class KMatrix: public Interfaced {

public:
  
  /**
   * Enum for the possible channels
   */
  enum Channels { PiPi, KPi, KEta, KEtaPrime};

public:

  /**
   * The default constructor.
   */
  KMatrix(FlavourInfo flavour=FlavourInfo(),
	  vector<Channels> channels=vector<Channels>(),
	  vector<Energy> poles=vector<Energy>());

  /**
   *   The quantum numbers of the K-matrix
   */
  FlavourInfo flavourInfo() const {
    return flavour_;
  };

  /**
   *  Compute the K-matrix for a given scale
   */
  virtual double K(Energy2 s) =0;
  
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KMatrix & operator=(const KMatrix &);

private:

  /**
   *   The quantum numbers for the K-matrix
   */
  FlavourInfo flavour_;
  
  /**
   *   The mesons in the various channels
   */
  vector<Channels> channels_;

  /**
   *  The positions of the poles
   */
  vector<Energy> poles_;
};

}

#endif /* Herwig_KMatrix_H */
