// -*- C++ -*-
//
// MamboDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MamboDecayer_H
#define HERWIG_MamboDecayer_H
//
// This is the declaration of the MamboDecayer class.
//

#include "HwDecayerBase.h"
#include "ThePEG/PDT/DecayMode.h"

namespace Herwig {
using namespace ThePEG;
  
  /**
   * The MamboDecayer class inherits from the Decayer class in 
   * ThePEG and implements the algorithm of R.Kleiss and 
   * W.J.Stirling NPB 385 (1992) 413-432 for massive multi-particle phase-space
   * decays 
   */
class MamboDecayer: public HwDecayerBase {

public:

  /**
   * The default constructor.
   */
  MamboDecayer() : _maxweight(10.), _a0(10,0.), _a1(10,0.) {}

  /**
   * Check if this decayer can perfom the decay for a particular mode
   * @param parent The decaying particle
   * @param children The decay products
   * @return true If this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   *  Perform the decay of the particle to the specified decay products
   * @param parent The decaying particle
   * @param children The decay products
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;


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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MamboDecayer> initMamboDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MamboDecayer & operator=(const MamboDecayer &);

private:

   /**
     *Set array of mometum to particles
     *@param mom   Momentum set to be distributed over phase-space
     *@param comEn The mass of the decaying particle
     *@return The weight of the configuration
     **/
  double calculateMomentum(vector<Lorentz5Momentum> & mom,
			   Energy comEn) const;

  /**
   * Set up the colour connections for the decay
   * @param parent The incoming particle
   * @param out The decay products
   */
  void colourConnections(const Particle & parent, 
			 ParticleVector & out) const;  

  /** @name Bessel Functions.*/
  //@{
  /**
   * Compute the values \f$K_0(x)/K_1(x)\f$ and it's derivative using
   * asymptotic expansion for large x values.
   * @param x The argument
   * @param f The value of the ratio
   * @param fp The value of the derivative ratio
   */
  void BesselFns(const long double x,
		 long double & f, long double & fp) const {
    assert(x>=0.);
    if( x < 10. ) {
      f = BesselK0(x)/BesselK1(x);
      fp = ( sqr(f)*x + f - x )/x;
    }
    else
      BesselIExpand(-x, f, fp);
  } 
  
  /**
   * Compute the values \f$I_0(x)/I_1(x)\f$ and it's derivative using
   * asymptotic expansion.
   * @param x The argument
   * @param f The value of the ratio
   * @param fp The value of the derivative ratio
   */
  void BesselIExpand(const long double x,
		     long double & f, long double & fp) const {
    long double y = 1./x;
    f = 1.+ y*(_a0[0] + y*(_a0[1] + y*(_a0[2] + y*(_a0[3] 
        + y*(_a0[4] + y*(_a0[5] + y*(_a0[6] + y*(_a0[7] 
        + y*(_a0[8] + y*_a0[9] )))))))));
    fp = -y*y*(_a1[0] + y*(_a1[1] + y*(_a1[2] + y*(_a1[3] 
        + y*(_a1[4] + y*(_a1[5] + y*(_a1[6] + y*(_a1[7] 
        + y*(_a1[8] + y*_a1[9] )))))))));
  }

  /**
   * Modified Bessel function of first kind \f$I_0(x)\f$.
   *@param x Argument of Bessel Function 
   **/
  long double BesselI0(const long double x) const {
    long double y,ans;
    if(x < 3.75) {
      y = sqr(x/3.75);
      ans = 1. + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 
          + y*(0.2659732 + y*(0.0360768+y*0.0045813)))));
    }
    else {
      y = (3.75/x);
      ans = (exp(x)/sqrt(x))*(0.39894228 + y*(0.01328592 
          + y*(0.00225319 + y*(-0.00157565 + y*(0.00916281 
          + y*(-0.02057706+y*(0.02635537+y*(-0.01647633+y*0.00392377))))))));
    }
    return ans;
  }
  
  /**
   *  Modified Bessel function of first kind \f$I_1(x)\f$.
   *@param x Argument of Bessel Function 
   **/
  long double BesselI1(const long double x) const {
    long double y,ans;
    if(x < 3.75) {
      y = sqr(x/3.75);
      ans = x*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934 
          + y*(0.02658733 + y*(0.00301532 + y*0.00032411))))));
    }
    else {
      y = 3.75/x;
      ans = (0.39894228 + y*(-0.03988024 + y*(-0.00362018 
          + y*(0.00163801 + y*(-0.01031555 + y*(0.02282967 
	  + y*(-0.02895312 + y*(0.01787654-y*0.00420059))))))))*(exp(x)/sqrt(x));
    }
    return ans;
  }
  
  /**
   * Modified Bessel function of second kind \f$K_0(x)\f$.
   * @param x Argument of Bessel Function 
   **/
  long double BesselK0(const long double x) const {
    long double y,ans;
    if(x <= 2.0) {
      y = x*x/4.0;
      ans = -log(x/2.0)*BesselI0(x) - 0.57721566 
          + y*(0.42278420 + y*(0.23069756 
          + y*(0.03488590 + y*(0.00262698 + y*(0.00010750+y*0.00000740)))));
    }
    else {
      y = 2.0/x;
      ans = (1.25331414 + y*(-0.07832358 + y*(+0.02189568 
          + y*(-0.01062446 + y*(0.00587872 
          + y*(-0.00251540 + y*0.00053208))))))*(exp(-x)/sqrt(x));
    }
    return ans;
  }
  
  /**
   * Modified Bessel function of second kind \f$K_1(x)\f$.
   * @param x Argument of Bessel Function 
   **/
  long double BesselK1(const long double x) const  {
    long double y,ans;
    if(x <= 2.0) {
      y = x*x/4.;
      ans = log(x/2.)*BesselI1(x) + (1./x)*(1. + y*(0.15443144 
          + y*(-0.67278579 + y*(-0.18156897 
          + y*(-0.01919402+y*(-0.00110404-(y*0.00004686)))))));
    }
    else {
      y = 2./x;
      ans = (exp(-x)/sqrt(x))*(1.25331414 + y*(0.23498619 
          + y*(-0.03655620 + y*(0.01504268 + y*(-0.00780353 
          + y*(0.00325614+y*(-0.00068245)))))));
    }
    return ans;
  }
  //@}

private:
  
  /**
   * Maximum weight
   */
  double _maxweight;

  /**
   * Store coefficents for aysymptotic expansion of \f$\frac{I_0}{I_1}\f$
   */
  vector<double> _a0;

  /**
   * Store data for aysymptotic expansion of the first derivative
   * \f$\frac{I_0}{I_1}\f$.
   */
  vector<double> _a1;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MamboDecayer. */
template <>
struct BaseClassTrait<Herwig::MamboDecayer,1> {
  /** Typedef of the first base class of MamboDecayer. */
  typedef Herwig::HwDecayerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MamboDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MamboDecayer>
  : public ClassTraitsBase<Herwig::MamboDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MamboDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MamboDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwMamboDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MamboDecayer_H */
