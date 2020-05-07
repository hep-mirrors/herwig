// -*- C++ -*-
#ifndef Herwig_MEDM2Mesons_H
#define Herwig_MEDM2Mesons_H
//
// This is the declaration of the MEDM2Mesons class.
//

#include "Herwig/MatrixElement/MEMultiChannel.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEDM2Mesons class implements dark matter annhilation via a vector current to
 * mesons at low energies
 *
 * @see \ref MEDM2MesonsInterfaces "The interfaces"
 * defined for MEDM2Mesons.
 */
class MEDM2Mesons: public MEMultiChannel {

public:

  /**
   * The default constructor.
   */
  MEDM2Mesons();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;
  //@}

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);

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
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   */
  virtual double me2(const int ichan) const;
  
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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEDM2Mesons & operator=(const MEDM2Mesons &) = delete;

private :

  /**
   *    Hadronic current etc
   */
  //@{
  /**
   * the hadronic current
   */
  WeakCurrentPtr current_;

  /**
   *  The matrix element
   */
  mutable ProductionMatrixElement me_;

  /**
   *  Map for the modes
   */
  map<int,int>  modeMap_;
  //@}

  /**
   *  DM
   */
  //@{
  /**
   *   Incoming Particles
   */
  PDPtr incomingA_, incomingB_;
  //@}

  /**
   * DM coupling to the dark mediator
   */
  Complex cDMmed_;

  /**                                                                                                                                                       
   * SM couplings to the dark mediator
   */
  vector<Complex> cSMmed_;

  /**
   * DM vector mediator
   */
  PDPtr Mediator_;
};

}

#endif /* Herwig_MEDM2Mesons_H */
