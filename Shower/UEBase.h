// -*- C++ -*-
#ifndef HERWIG_UEBase_H
#define HERWIG_UEBase_H
//
// This is the declaration of the UEBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/StandardXComb.fh"
#include "UEBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Abstract base class used to minimize the dependence between
 * MPIHandler and all Shower classes.
 *
 * \author Manuel B\"ahr
 *
 * @see \ref UEBaseInterfaces "The interfaces"
 * defined for UEBase.
 */
class UEBase: public Interfaced {

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /** @name Virtual functions used for the generation of additional
      interactions . */
  //@{
  /**
   * Some initialization code eventually.
   */
  virtual void initialize() {}

  /**
   * Return true or false depending on the generator setup. 
   */
  virtual bool beamOK() const = 0;

  /**
   * Return true or false depending on whether soft interactions are enabled.
   */
  virtual bool softInt() const {return false;}

  /**
   * Return the value of the pt cutoff.
   */
  virtual Energy Ptmin() const = 0;

  /**
   * Return the slope of the soft pt spectrum. Only necessary when the
   * soft part is modelled.
   */
  virtual InvEnergy2 beta() const {return ZERO;}

  /**
   * Some finalize code eventually.
   */
  virtual void finalize() {}

  /**
   * Clean up method called after each event.
   */
  virtual void clean() {}

  /**
   * Return the number of different hard processes. Use 0 as default to
   * not require implementation.
   */
  virtual unsigned int additionalHardProcs() const {return 0;}

  /**
   * return the hard multiplicity of process i. Can't be constant in my
   * case because drawing from the probability distribution also
   * specifies the soft multiplicity that has to be stored....
   */
  virtual unsigned int multiplicity(unsigned int i=0) = 0;

  /**
   * Generate a additional interaction for ProcessHandler sel. Method
   * can't be const because it saves the state of the underlying XComb
   * object on it's way.
   */
  virtual tStdXCombPtr generate(unsigned int sel=0) = 0;

  /**
   * Return the type of algorithm. 
   */
  virtual int Algorithm() const = 0;

  /**
   * Return the value of the hard Process pt cutoff for vetoing.
   */
  virtual Energy PtForVeto() const = 0;

  /**
   * Return the fraction of colour disrupted subprocesses. Use default 0
   * so that it is not required to implement.
   */
  virtual double colourDisrupt() const {return 0.0;}

  /**
   * Return the soft multiplicity. Use 0 as default to not require
   * implementation.
   */
  virtual unsigned int softMultiplicity() const {return 0;} 
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEBase & operator=(const UEBase &);

};

}

#endif /* HERWIG_UEBase_H */
