// -*- C++ -*-
#ifndef Herwig_GGPsiSplitFn_H
#define Herwig_GGPsiSplitFn_H
//
// This is the declaration of the GGPsiSplitFn class.
//

#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The documentation of the GGPsiSplitFn class implements the colour singlet spliting function for \f$g\to g J\Psi,\Upsilon\f$
 *
 * @see \ref GGPsiSplitFnInterfaces "The interfaces"
 * defined for GGPsiSplitFn.
 */
class GGPsiSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   * The default constructor.
   */
  GGPsiSplitFn();
  
public:
  
  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    for(unsigned int ix=0;ix<2;++ix) {
      if(ids[ix]->id()!=ParticleID::g) return false;
    }
    return true;
  }

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * Value of the energy fraction and value of the scale for the veto algorithm
   * @param iopt The option for calculating z
   * @param ids The PDG codes of the particles in the splitting
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   * @param t1 The starting valoe of the scale
   * @param enhance The radiation enhancement factor
   * @param identical Whether or not the outgoing particles are identical
   * @param t_main rerurns the value of the energy fraction for the veto algorithm
   * @param z_main returns the value of the scale for the veto algorithm
   */
  virtual void guesstz(Energy2 t1,unsigned int iopt, const IdList &ids,
		       double enhance,bool ident,
		       double detune, Energy2 &t_main, double &z_main);
  
  /**
   * The concrete implementation of the splitting function, \f$P(z,t)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   * @param rho The spin density matrix
   */
  double P(const double z, const Energy2 t,
	   const IdList & , const bool, const RhoDMatrix &) const {
    double r = 4.*sqr(m_)/t;
    return  I(z,r)/z/(1.-z);
  }

  /**
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  double overestimateP(const double z, const IdList &) const {
    return maxP_/z/(1.-z); 
  }

  /**
   * The concrete implementation of the
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,t)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   * @param rho The spin density matrix
   */
  double ratioP(const double z, const Energy2 t,
		const IdList & , const bool, const RhoDMatrix &) const {
    double r = 4.*sqr(m_)/t;
    double ratio  = I(z,r)/maxP_;
    if(ratio>1.) cerr << "Overestimate violated " <<ratio << " " << z << " " << t/GeV2 << " " << y_ << "\n";
    return ratio;
  }
  
  /**
   * The concrete implementation of the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  double integOverP(const double z, const IdList & ,
		    unsigned int PDFfactor=0) const {
    assert(PDFfactor==0 && z>0. && z<1.);
    return maxP_*log(z/(1.-z)); 
  }

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  double invIntegOverP(const double r, const IdList & ,
			       unsigned int PDFfactor=0) const {
    assert(PDFfactor==0);
    return 1./(1.+exp(-r/maxP_));
  }

  /**
   * Method to calculate the azimuthal angle for forward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  vector<pair<int,Complex> >
  generatePhiForward(const double z, const Energy2 t, const IdList &,
		     const RhoDMatrix & rho) {
    assert(rho.iSpin()==PDT::Spin1);
    double r = sqr(m_)/(t-4.*sqr(m_));
    double on  = 0.5*(1.-2.*z*(1.-z)) - 4.*z*(1.-2.*z)*r + 16.*sqr(z*r);
    double off = z*(1.-z) + 4.*z*(1.-2.*z)*r - 16.*sqr(z*r);
    double max = on+ 2.*abs(rho(0,2))*off;
    vector<pair<int, Complex> > output;
    output.reserve(3);
    output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*on/max));
    output.push_back(make_pair(-2, rho(0,2)*off/max));
    output.push_back(make_pair( 2, rho(2,0)*off/max));
    return output;
  }

  /**
   * Method to calculate the azimuthal angle for backward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  vector<pair<int,Complex> >
  generatePhiBackward(const double z, const Energy2, const IdList &,
		      const RhoDMatrix & rho) {
    assert(false);
    assert(rho.iSpin()==PDT::Spin1);
    double diag = sqr(1 - (1 - z)*z)/(1 - z)/z;
    double off  = (1.-z)/z;
    double max  = 2.*abs(rho(0,2))*off+diag;
    vector<pair<int, Complex> > output;
    output.reserve(3);
    output.push_back(make_pair( 0, (rho(0,0)+rho(2,2))*diag/max));
    output.push_back(make_pair( 2,-rho(0,2)           * off/max));
    output.push_back(make_pair(-2,-rho(2,0)           * off/max));
    return output;
  }
  
  
  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   */
  DecayMEPtr matrixElement(const double z, const Energy2 t, 
			   const IdList &, const double phi, bool) {
    // calculate the kernal
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
    Complex ii(0.,1.);
    double r  = sqr(m_)/(t-4.*sqr(m_));
    Complex phase = exp(2.*ii*phi);
    (*kernal)(0,0,0) =  z*(1.+4.*r);
    (*kernal)(2,2,0) = -(*kernal)(0,0,0) ;
    (*kernal)(0,2,0) = -(1.-z-4.*z*r)/phase;
    (*kernal)(2,0,0) = -conj((*kernal)(0,2,0));
    return kernal;
  }
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

  double alphaSVetoRatio(Energy2 pt2, double factor) const {
    factor *= ShowerHandler::currentHandler()->renormalizationScaleFactor();
    return pow(alpha()->ratio(pt2, factor),3);
  }

  /**
   *  Phase Space veto member to implement the \f$\Theta\f$ function as a veto
   *  so that the emission is within the allowed phase space.
   * @param t  The scale
   * @param maxQ2 The maximum virtuality
   * @return true if vetoed
   */
  virtual bool PSVeto(const Energy2 t) {
    double r = 4.*sqr(m_)/z()/(1.-z())/t;
    if ( y_>0.5*(1.+r) || y_<0.5*(r+sqr(1.-z()))/(1.-z()) ) return true;
    return Sudakov1to2FormFactor::PSVeto(t);
  }
  
protected:

  /**
   *  Integrand for the splitting function
   */
  double I(double zz, double r) const {
    double z = 1.-zz;
    double r2(sqr(r)),r3(r*r2);
    double y2(sqr(y_)),y3(y2*y_),y4(y3*y_),y5(y4*y_),y6(y5*y_);
    double f0 = r2*(1 + r)*(3 + 12*r + 13*r2) - 16*r2*(1 + r)*(1 + 3*r)*y_
      - 2*r*(3 - 9*r - 21*r2 + 7*r3)* y2 + 8*r*(4 + 3*r + 3*r2)*y3
      - 4*r*(9 - 3*r - 4*r2)*y4 - 16*(1 + 3*r + 3*r2)*y5 + 8*(6 + 7*r)*y6 - 32*y_*y6;
    double f1 = -2*r*(1 + 5*r + 19*r2 + 7*r3)*y_ + 96*r2*(1 + r)*y2
      + 8*(1 - 5*r - 22*r2 - 2*r3)*y3 + 16*r*(7 + 3*r)*y4 - 8*(5 + 7*r)*y5 + 32*y6;
    double f2 = r*(1 + 5*r + 19*r2 + 7*r3) - 48*r2*(1 + r)*y_ - 4*(1 - 5*r - 22*r2 - 2*r3)*y2
      - 8*r*(7 + 3*r)*y3 + 4*(5 + 7*r)*y4 - 16*y5;
    double g0 = (1 - r)*r3*(3 + 24*r + 13*r2) - 4*r3*(7 - 3*r - 12*r2)*y_
      - 2*r3*(17 + 22*r - 7*r2)*y2 + 4*r2*(13 + 5*r - 6*r2)*y3
      - 8*r*(1 + 2*r + 5*r2 + 2*r3)*y4 - 8*r*(3 - 11*r - 6*r2)*y5 + 8*(1 - 2*r - 5*r2)*y6;
    double g1 = -2*(1 - r)*r2*(1 + r)*(1 + 7*r)*y_ + 8*(1 - 4*r)*r2*(1 + 3*r)*y2
      + 4*r*(1 + 10*r + 57*r2 + 4*r3)* y3 - 8*r*(1 + 29*r + 6*r2)*y4 - 8*(1 - 8*r - 5*r2)*y5;
    double g2 = (1 - r)*r2*(1 + r)*(1 + 7*r) - 4*(1 - 4*r)*r2*(1 + 3*r)*y_
      - 2*r*(1 + 10*r + 57*r2 + 4*r3)* y2 + 4*r*(1 + 29*r + 6*r2)*y3 + 4*(1 - 8*r - 5*r2)*y4;
    double c1 = f0+z*f1+sqr(z)*f2;
    double c2 = g0+z*g1+sqr(z)*g2;
    double root = sqrt(-r + y2);
    return (c1 + ((1 + r - 2*y_)*c2*log((-r + y_ + root)/(-r + y_ - root)))/(2.*(-r + y_)*root))
      /(sqr(1 - y_)*sqr(-r + y_)*sqr(-r + y2));
  }

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
  GGPsiSplitFn & operator=(const GGPsiSplitFn &) = delete;

private:

  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy3 O1_;

  /**
   *  The quark mass
   */
  Energy m_;

  /**
   *  Option for the quark mass
   */
  unsigned int massOpt_;

  /**
   *  Maximum value for the overestimate
   */
  double maxP_;

  /**
   *  Auxillary y variable
   */
  double y_;

};

}

#endif /* Herwig_GGPsiSplitFn_H */
