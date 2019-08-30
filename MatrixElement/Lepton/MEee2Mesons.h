// -*- C++ -*-
#ifndef Herwig_MEee2Mesons_H
#define Herwig_MEee2Mesons_H
//
// This is the declaration of the MEee2Mesons class.
//

#include "Herwig/MatrixElement/MEMultiChannel.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2Mesons class implements \f$e^+e^-\f$ annhilation to mesons at low energy
 * using hadronic currents.
 *
 * @see \ref MEee2MesonsInterfaces "The interfaces"
 * defined for MEee2Mesons.
 */
class MEee2Mesons: public MEMultiChannel {

public:

  /**
   * The default constructor.
   */
  MEee2Mesons();

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
  
  /**
   * Return the matrix element squared for a given mode and phase-space channel,
   * with the helicities amplitudes
   * @param ichan The channel we are calculating the matrix element for. 
   */
  virtual double helicityME(const int ichan, const cPDVector & particles,
			    const vector<Lorentz5Momentum> & momenta) const;

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

  /**
   *  Set up the flavours allowed in the current
   */
  void setFlavour();
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2Mesons & operator=(const MEee2Mesons &) = delete;

private :

  /**
   *  Option for the flavour of the particles in the current
   */
  unsigned int flavOpt_;
  
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

  /**
   *  The flavour of the hadronic system
   */
  FlavourInfo flavour_;
};

}

#endif /* Herwig_MEee2Mesons_H */
