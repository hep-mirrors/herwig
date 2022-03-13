// -*- C++ -*-
#ifndef Herwig_GtoG3P1SplitFn_H
#define Herwig_GtoG3P1SplitFn_H
//
// This is the declaration of the GtoG3P1SplitFn class.
//

#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/MatrixElement/Onium/OniumParameters.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GtoG3P1SplitFn class implements the colour singlet splitting function for \f$g\to g \chi_{0(c,b)}\f$.
 *
 * @see \ref GtoG3P1SplitFnInterfaces "The interfaces"
 * defined for GtoG3P1SplitFn.
 */
class GtoG3P1SplitFn: public Sudakov1to2FormFactor {

public:
  
  /**
   * The default constructor.
   */
  GtoG3P1SplitFn() : O1_(0.573*GeV*GeV2*GeV2), state_(ccbar), n_(1), m_(1.2*GeV), fixedAlphaS_(-1.)
  {}
    
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
    // check onium state
    int iq=4+state_;
    long idtest = iq*110+20003 + (n_-1)*100000;
    if(ids[2]->id() != idtest) return false;
    // charge conservation
    if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge()) return false;
    // looks OK
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
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  double overestimateP(const double z, const IdList &) const {
    return 1/z + 1/(1.-z); 
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
    double r = sqr(m_)/(t-4.*sqr(m_));
    double ratio = 1.  - 2.*(1.-z)*z + 8*r*(1 + z*(-1 + 2*z))  + 16*sqr(r)*(1 + 2*z)  + 128*pow(r,3)*(1 - 2*z)*z - 512*pow(r,4)*sqr(z);
    ratio /= pOver_;
    if(ratio>1.) cerr << "problem in GtoG3P1 ratio violated " << ratio << " " << pOver_ << "\n";
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
    return log(z/(1.-z)); 
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
    return 1./(1.+exp(-r));
  }
  //@}

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
    double r  = sqr(m_)/(t-4.*sqr(m_));
    double on  = 1.-2.*(1.-4.*r)*z*(1.-z-4.*r*z);
    double off = 2.*(1.-4.*r)*z*(1.-z-4.*r*z);
    double max = on+ 2.*abs(rho(0,2)*off);
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
  generatePhiBackward(const double, const Energy2, const IdList &,
		      const RhoDMatrix & rho) {
    assert(false);
    assert(rho.iSpin()==PDT::Spin1);
    vector<pair<int, Complex> > output;
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
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
    Complex ii(0.,1.);
    double r  = sqr(m_)/(t-4.*sqr(m_));
    double rt = sqrt(r*z*(1.-z-4*r*z));
    Complex phase = exp(ii*phi);
    (*kernal)(0,0,0) = (2*sqrt(2)*phase*rt)/(-1 + z);
    (*kernal)(0,0,1) = (z*(-1 + z + 4*r*(1 + z)))/(-1 + z);
    (*kernal)(0,0,2) = (-2*sqrt(2)*z*rt)/(phase*(-1 + z));
    (*kernal)(0,2,0) = (2*sqrt(2)*rt)/phase;
    (*kernal)(0,2,1) = (-1 + z + 4*r*z)/sqr(phase);
    (*kernal)(0,2,2) = 0.;
    (*kernal)(2,0,0) = 0.;
    (*kernal)(2,0,1) = -(sqr(phase)*(-1 + z + 4*r*z));
    (*kernal)(2,0,2) = 2*sqrt(2)*phase*rt;
    (*kernal)(2,2,0) = (-2*sqrt(2)*phase*z*rt)/(-1 + z);
    (*kernal)(2,2,1) = -((z*(-1 + z + 4*r*(1 + z)))/(-1 + z));
    (*kernal)(2,2,2) = (2*sqrt(2)*rt)/(phase*(-1 + z));
    return kernal;
  }

protected:

  double alphaSVetoRatio(Energy2 pt2, double factor) const {
    if(fixedAlphaS_<0.) {
      factor *= ShowerHandler::currentHandler()->renormalizationScaleFactor();
      return sqr(alpha()->ratio(pt2, factor));
    }
    else {
      return 1.;
    }
  }

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

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
  GtoG3P1SplitFn & operator=(const GtoG3P1SplitFn &) = delete;

private:
  
  /**
   *  Access to the parameters for the quarkonium states
   */
  OniumParametersPtr params_;

  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy5 O1_;

  /**
   *  Type of state
   */
  OniumState state_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;

  /**
   *  The quark mass
   */
  Energy m_;

  /**
   *  Fixed value of \f$\alpha_S\f$
   */
  double fixedAlphaS_;

  /**
   *  Max of the spltting function
   */
  static const double pOver_;
  
};

}

#endif /* Herwig_GtoG3P1SplitFn_H */
