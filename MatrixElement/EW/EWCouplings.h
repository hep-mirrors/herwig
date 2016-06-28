// -*- C++ -*-
#ifndef Herwig_EWCouplings_H
#define Herwig_EWCouplings_H
//
// This is the declaration of the EWCouplings class.
//

#include "ThePEG/Interface/Interfaced.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "EWCouplings.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the EWCouplings class.
 *
 * @see \ref EWCouplingsInterfaces "The interfaces"
 * defined for EWCouplings.
 */
class EWCouplings: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EWCouplings(unsigned int loops=2, unsigned int steps=200,
	      Energy highScale=10.*TeV, Energy lowScale=10.*GeV);

  /**
   * The destructor.
   */
  virtual ~EWCouplings();
  //@}

  /**
   *  Initialize the running couplings
   */
  void initialize();

  /**
   * Number of dynamical quarks at $\log\mu = x$ (in GeV)
   * N.B.Integrate out top quark at Mz, not at Mt.
   */
  unsigned int numberOfFlavours(double x) {
    return x >= std::log(ewScale_/GeV) ? 6 : 5;
  }
   
public:
   
  /**
   *  Set whether or not to include \f$SU(3)\f$ running
   */
  void SU3(bool includeSU3) {includeSU3_ = includeSU3;}

  /**
   *  Whether or not to include \f$SU(3)\f$ running
   */
  bool SU3() { return includeSU3_;}
    
  /**
   *  Set whether or not to include EW running
   */  
  void EW(bool includeEW) {includeEW_ = includeEW;}

  /**
   *  Whether or not to include EW running
   */
  bool EW() { return includeEW_;}
   
  /**
   *  alpha for the U1 gauge coupling at energy mu (in GeV):
   */
  double a1(Energy mu) {
    if (includeEW_) {
      if (mu>=ewScale_) {
	return (3.0/5.0)*interpolate(log(mu/GeV),1);
      }
      return interpolateLow(log(mu/GeV),1);
    }
    else
      return 0.0;
  }

  /**
   *  alpha for the SU2 gauge coupling at energy mu (in GeV):
   */
  double a2(Energy mu) {
    if (includeEW_) {
      if (mu<ewScale_) {
	return 0.0;
      }
      else
	return interpolate(log(mu/GeV),2);
    }
    else
      return 0.0;
  }

  /**
   *  alpha for the SU3 gauge coupling at energy mu (in GeV):
   */
  double a3(Energy mu) {
    if(includeSU3_) {
      if (mu>=ewScale_) {
	return interpolate(log(mu/GeV),3);
      }
      else {
	return interpolateLow(log(mu/GeV),2);
      }
    }
    else
      return 0.0;
  }

  /**
   *  alpha for EM
   */
  double aEM(Energy mu) {
    if (includeEW_) {
      if (mu<=ewScale_) {
	return interpolateLow(log(mu/GeV),1);
      }
      else {
	double alpha1=a1(mu);
	double alpha2=a2(mu);
	return alpha1*alpha2/(alpha1+alpha2);
      }
    }
    return 0.0;
  }

  double aS(Energy mu) {
    if(includeSU3_) {
      if (mu<=ewScale_) {
	return interpolateLow(log(mu/GeV),2);
      }
      else {
	return interpolate(log(mu/GeV),3);
      }
    }
    else return 0.0;
  }

  /**
   *  Top quark Yukawa coupling
   */
  double y_t(Energy mu) {
    if (includeEW_) {
      if(mu<ewScale_)
	return 0.0;
      else
	return interpolate(log(mu/GeV),4);
    }
    else
      return 0.0;
  }

  /**
   * Quartic scalar coupling lambda (normalization different, hence factor of 2):
   */
  double lambda(Energy mu) {return 2.0*interpolate(log(mu/GeV),5);}
  
  /**
   *  VEV
   */
  Energy vev(Energy mu) {return interpolate(log(mu/GeV),6)*GeV;}
      
  /**
   * \f$\lambda_t\f$
   */
  double lambda_t(Energy mu) {return y_t(mu);}
      
  /**
   *  Z couplings
   */
  //@{
  double g_Lu(Energy mu) {return 0.5-(2.0/3.0)*Sin2thW(mu);}
  double g_Ld(Energy mu) {return -0.5+(1.0/3.0)*Sin2thW(mu);}
  double g_Le(Energy mu) {return (-1.0/2.0)+Sin2thW(mu);}
  double g_Lnu(Energy ) {return (0.5);}
  double g_Ru(Energy mu) {return (-2.0/3.0)*Sin2thW(mu);}
  double g_Rd(Energy mu) {return (1.0/3.0)*Sin2thW(mu);}
  double g_Re(Energy mu) {return Sin2thW(mu);}
  double g_W(Energy mu) {return Cos2thW(mu);}
  double g_phiPlus(Energy mu) {return 0.5-Sin2thW(mu);}
  //@}

  /**
   *  \f\cos^2\theta_W\f$
   */
  double Cos2thW(Energy ) {
    //\todo why this value?
    //return mW()*mW()/(mZ()*mZ());
    return (1.-0.2314);
  }

  /**
   *  \f\sin^2\theta_W\f$
   */
  double Sin2thW(Energy mu) {return 1.-Cos2thW(mu);}

  double aW(Energy mu) {return a2(mu);}

  double aZ(Energy mu) {
    //return a2(mu)/Cos2thW(mu); // Same thing, actually
    return a1(mu)+a2(mu);
  }

public: 

  /**
   *  Masses of the Standard Model particles
   */
  //@{
  /**
   *  Z mass
   */
  Energy mZ() const { return mZ_;}
 
  /**
   *  W Mass
   */
  Energy mW() const { return mW_;}

  /**
   *  Top quark mass
   */
  Energy mT() const {return mT_;}

  Energy mTatmZ() { return mT_;}

  /**
   *  Higg boson mass
   */
  Energy mH() const {return mH_;}
  //@}

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
   *  Set the initial value of the couplings
   */
  void initializeCouplings(vector<Complex> & y);

  /**
   * Assigns numerical values to beta functions 
   * Takes in a point x = log(mu) and the values of y[i] at x and assigns dydx[i] the 
   * value beta[i](mu).  The function Derivs farms out the plugging in to the three
   * functions BetaGauge, BetaYukawa, and BetaHiggs, which evaluates the beta functions
   * for the gauge couplings, yukawa matrices, and higgs quartic coupling/vev, respectively.
   */
  void derivatives(const double x, vector< Complex > &y,
		   vector< Complex > & dydx);

  /**
   * Beta function for the gauge interactions
   */
  void betaGauge(const double x, vector<Complex> &y, vector<Complex> &dydx);

  /**
   * Beta function for the gauge interactions at low scales
   */
  void lowBetaGauge(const double x, vector<Complex> &y, vector<Complex> &dydx);

  /**
   * Beta function for the Yukawa interactions
   */
  void betaYukawa(const double x, vector<Complex> &y, vector<Complex> &dydx);

  /**
   * Beta function for the Higgs interactions
   */
  void betaHiggs(const double x, vector<Complex> &y, vector<Complex> &dydx);

  /**
   *  Update the couplings using 4-th order RK
   * Takes in a vector y[i] of function values at a point x and the values of the 
   * first derivatives dydx[i] ( = dy[i]/dx ) alon with a step size stepsize. The 
   * function then defines assigns the value of y[i](x+stepsize) to the variable yout[i].
   * (Adapted from sample code in Numerical Recipes in C++ Press, Teukolsky, et. al.)
   */
  void RK4(vector<Complex> & y, vector<Complex> &dydx, const double x, 
	   const double stepsize, vector<Complex> &yout);

  /**
   *  Initialize the low energy parameters
   */
  void initializeLow();

  /**
   *  Interpolate the table,  t = ln(mu)
   */
  double interpolate(double t, int paramIndex);

  /**
   *  Interpolate the tabel, t = ln(mu)
   */
  double interpolateLow(double t, int paramIndex);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EWCouplings & operator=(const EWCouplings &);

private:

  /**
   *  Electoweak Scale
   */
  Energy ewScale_;

  /**
   *  High Scale
   */
  Energy highScale_;

  /**
   *  Low Scale
   */
  Energy lowScale_;

  /**
   *  Whether or not to include SU3
   */
  bool includeSU3_;

  /**
   *  Whether or not to include EW
   */
  bool includeEW_;

  /**
   *  Whether or not the running couplings have been initialized
   */
  bool initialized_;

  /**
   *  Masses of Standard Model particles
   */
  //@{
  /**
   *  Mass Choice
   */
  bool massChoice_;

  /**
   *   Z mass 
   */
  Energy mZ_;

  /**
   *   W mass 
   */
  Energy mW_;

  /**
   *  Top mass
   */
  Energy mT_;

  /**
   *  Higgs boson mass
   */
  Energy mH_;
  //@}

  /**
   *  Number of loops
   */
  unsigned int loops_;

  /**
   *  Number of steps for Runga-Kutta (High scale)
   */
  unsigned int highSteps_;

  /**
   *  Number of steps for Runga-Kutta (Low scale)
   */
  unsigned int lowSteps_;

  /**
   *  Matrix to store the parameters
   */
  boost::numeric::ublas::matrix<double> table_;

  /**
   *  Matrix to store the low energy parameters.
   *  This will hold only aEM and a3 at mu<ewScale
   */
  boost::numeric::ublas::matrix<double> lowTable_;

};

}

#endif /* Herwig_EWCouplings_H */
