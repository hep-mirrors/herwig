// -*- C++ -*-
#ifndef Herwig_SoftSudakov_H
#define Herwig_SoftSudakov_H
//
// This is the declaration of the SoftSudakov class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Utilities/GSLIntegrator.h"
// work around a Boost 1.64 bug where ublas headers would fail otherwise
#include <boost/version.hpp>
#if (BOOST_VERSION / 100 >= 1064)
#include <boost/serialization/array_wrapper.hpp>
#endif
#include <boost/numeric/ublas/matrix.hpp>
#include "EWProcess.h"
#include "SoftSudakov.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SoftSudakov class.
 *
 * @see \ref SoftSudakovInterfaces "The interfaces"
 * defined for SoftSudakov.
 */
class SoftSudakov: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SoftSudakov();

  /**
   * The destructor.
   */
  virtual ~SoftSudakov();
  //@}

public:

  /**
   *  Low energy soft evolution
   */
  boost::numeric::ublas::matrix<Complex>
  lowEnergyRunning(Energy EWScale, Energy lowScale, 
		   Energy2 s, Energy2 t, Energy2 u,
		   Herwig::EWProcess::Process process,
		   unsigned int iswap);
  
  /**
   *  Evalaute the high energy running as a matrix
   */
  boost::numeric::ublas::matrix<Complex>
  highEnergyRunning(Energy highScale, Energy EWScale,
		    Energy2 s, Energy2 t, Energy2 u,
		    Herwig::EWProcess::Process process,
		    unsigned int iswap);

  /**
   *  Number of operators for the broken theory at low energy
   */
  unsigned int numberBrokenGauge(Herwig::EWProcess::Process process);

  /**
   *  Number of operators for the unbroken theory at high energy
   */
  unsigned int numberGauge(Herwig::EWProcess::Process process);

protected:
      
  /**
   *  Evaluate the soft evolution 
   */
  boost::numeric::ublas::matrix<Complex> evaluateSoft(boost::numeric::ublas::matrix<Complex> & G3,
						      boost::numeric::ublas::matrix<Complex> & G2,
						      boost::numeric::ublas::matrix<Complex> & G1,
						      Energy mu_h, Energy mu_l, bool high);

public:

  /**
   *  The operator to be integrated
   */
  InvEnergy operator ()(Energy mu) const;
  /** Argument type for GaussianIntegrator */
  typedef Energy ArgType;
  /** Return type for GaussianIntegrator */
  typedef InvEnergy ValType;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SoftSudakov & operator=(const SoftSudakov &);

private:

  /**
   *  Order for K
   */
  unsigned int K_ORDER_;

  /**
   *  Integrator
   */
  GSLIntegrator integrator_;

  /**
   *  Whether doing the high or low scale contribution
   */
  bool high_;

  /**
   *  Column
   */
  unsigned int row_; 

  /**
   *  Row
   */
  unsigned int col_; 

  /**
   *
   */
  bool real_; 

  /**
   *
   */
  boost::numeric::ublas::matrix<Complex> G1_;

  /**
   *
   */
  boost::numeric::ublas::matrix<Complex> G2_;

  /**
   *
   */
  boost::numeric::ublas::matrix<Complex> G3_;
};

}

#endif /* Herwig_SoftSudakov_H */
