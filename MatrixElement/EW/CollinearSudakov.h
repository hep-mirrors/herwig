// -*- C++ -*-
#ifndef Herwig_CollinearSudakov_H
#define Herwig_CollinearSudakov_H
//
// This is the declaration of the CollinearSudakov class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Utilities/GSLIntegrator.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "EWProcess.h"
#include "CollinearSudakov.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  Struct for the wavefunction corrections
 */
struct WaveFunctionCorrections {
  Complex RW;
  Complex RA;
  Complex RZ;
  Complex RAtoZ;
  Complex RZtoA;
  Complex RPhi;
  Complex EW;
  Complex EZ;
  Complex RPhi3;
  Complex RH;
  Complex tLuLDiff;
  Complex bLdLDiff;
  Complex tRuRDiff;
  
  // The following are constants from parameter integrals:
  
  Complex fFW0;
  Complex fF0W;
  Complex fFZZ;
  Complex aHH;
  Complex aZZ;
  Complex aW0;
  Complex a0W;
  Complex bHH;
  Complex bZZ;
  Complex cHH;
  Complex cZZ;
  Complex cW0;
  Complex atHH;
  Complex atZZ;
  Complex atW0;
  Complex at0W;
  Complex ctHH;
  Complex ctZZ;
  Complex ctW0;
  Complex btHH;
  Complex btZZ;
  
  Complex fs10;
  Complex fs1ZW;
  Complex fsWZWZ;
  Complex fsZW1;
  Complex fs01;
  Complex fsHW1;
  Complex fsHZ1;
  Complex fs1HW;
  Complex fs1HZ;
};

/**
 * Here is the documentation of the CollinearSudakov class.
 *
 * @see \ref CollinearSudakovInterfaces "The interfaces"
 * defined for CollinearSudakov.
 */
class CollinearSudakov: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CollinearSudakov();

  /**
   * The destructor.
   */
  virtual ~CollinearSudakov();
  //@}

public:
  
  /**
   *  Evalaute the electroweka matching as a matrix
   */
  boost::numeric::ublas::matrix<Complex>
  electroWeakMatching(Energy EWScale, Energy2 s,
		      Herwig::EWProcess::Process process,
		      bool oneLoop);

public:

  /**
   *  Make plots for tests
   */
  void makePlots();

protected:

  /**
   *   Evaluate the high scale contributions
   */
  void evaluateHighScale(Energy highScale, Energy EWScale, Energy2 S);

  /**
   *   Evaluate the low scale contributions
   */
  void evaluateLowScale(Energy EWScale, Energy lowScale, Energy2 S);

  /**
   *  Evaluate the matching
   */
  void evaluateMatching(Energy EWScale,Energy2 S);

public:

  /**
   *  The operator to be integrated
   */
  InvEnergy operator ()(Energy mu) const {
    if(high_) return highScaleIntegrand(mu);
    else      return  lowScaleIntegrand(mu);
  }
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

  /**
   *  The integral of the high scale part of the Sudakov
   */
  Complex highScaleIntegral(bool SU3, bool SU2, double Y,
			    Energy2 s, Energy mu_h, Energy mu_l, bool fermion, 
			    bool longitudinal, double yukFactor);

  /**
   *  the integral of the low scale part of the Sudakov
   */
  Complex  lowScaleIntegral(bool SU3, double Q, Energy2 s, 
			    Energy mu_h, Energy mu_l, bool fermion, 
			    double boostFactor);

protected:

  /**
   *   High-scale integrand
   */
  InvEnergy highScaleIntegrand(Energy mu) const;

  /**
   *   Low-scale integrand
   */
  InvEnergy  lowScaleIntegrand(Energy mu) const;

  /**
   *  Calculate the wavefunction corrections
   */
  WaveFunctionCorrections waveFunctionCorrections(Energy EWScale);

  /**
   *  Collinear matiching for W
   */
  Complex CollinearDw(Energy2 s, Energy EWScale);

  /**
   *  Collinear matching for Z
   */
  Complex CollinearDz(Energy2 s, Energy EWScale);

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
  CollinearSudakov & operator=(const CollinearSudakov &);
private:

  /**
   *  Parameters for the integrand
   */
  //@{
  /**
   *  Whether high or low scale
   */
  bool high_;

  /**
   *  Whether real of imaginary part
   */
  bool real_;

  /**
   * \f$SU(3)\f$
   */
  bool SU3_;

  /**
   *
   */
  bool SU2_;

  /**
   *
   */
  double Y_;

  /**
   *
   */
  Energy2 s_;

  /**
   *
   */
  bool fermion_; 

  /**
   *
   */
  bool longitudinal_;

  /**
   *
   */
  double yukFactor_;

  /**
   *
   */
  double boostFactor_;

  /**
   *
   */
  double Q_;
  //@}

  /**
   * Parameters
   */
  //@{
  /**
   *  Order for the K terms
   */
  int K_ORDER_;

  /**
   *  Order for the B terms
   */
  int B_ORDER_;
  //@}

  /**
   *  Integrator
   */
  GSLIntegrator integrator_;

private:

  /**
   *   Storage of the high scale pieces
   */
  //@{
  /**
   *
   */
  Complex highColW_;

  /**
   *
   */
  Complex highColB_;

  /**
   *
   */
  Complex highColG_;

  /**
   *
   */
  Complex highColQ_;

  /**
   *
   */
  Complex highColQt_;

  /**
   *
   */
  Complex highColU_;

  /**
   *
   */
  Complex highColtR_;

  /**
   *
   */
  Complex highColD_;

  /**
   *
   */
  Complex highColL_;

  /**
   *
   */
  Complex highColE_;

  /**
   *
   */
  Complex highColPhi_;
  //@}

  /**
   *  Storage of the low scale pieces
   */
  //@{
  /**
   *
   */
  complex<double> lowColW_;

  /**
   *
   */
  Complex lowColA_;

  /**
   *
   */
  Complex lowColG_;

  /**
   *
   */
  Complex lowColU_;

  /**
   *
   */
  Complex lowColt_;

  /**
   *
   */
  Complex lowColD_;

  /**
   *
   */
  Complex lowColE_;
  //@}

  /**
   *  Storage of the matching parameters
   */
  //@{
  /**
   *
   */
  Complex ULcollinearCorr_;

  /**
   *
   */
  Complex DLcollinearCorr_;

  /**
   *
   */
  Complex URcollinearCorr_;

  /**
   *
   */
  Complex DRcollinearCorr_;

  /**
   *
   */
  Complex tLcollinearCorr_;

  /**
   *
   */
  Complex tRcollinearCorr_;

  /**
   *
   */
  Complex bLcollinearCorr_;

  /**
   *
   */
  Complex nuLcollinearCorr_;

  /**
   *
   */
  Complex ELcollinearCorr_;

  /**
   *
   */
  Complex ERcollinearCorr_;

  /**
   *
   */  
  Complex WtoWcollinearCorr_;

  /**
   *
   */
  Complex WtoZcollinearCorr_;

  /**
   *
   */
  Complex WtoAcollinearCorr_;

  /**
   *
   */
  Complex BtoZcollinearCorr_;

  /**
   *
   */
  Complex BtoAcollinearCorr_;

  /**
   *
   */
  Complex PhitoWcollinearCorr_;

  /**
   *
   */
  Complex PhitoZcollinearCorr_;

  /**
   *
   */
  Complex PhitoHcollinearCorr_;

  /**
   *
   */
  Complex GcollinearCorr_;
  //@}
};

}

#endif /* Herwig_CollinearSudakov_H */
