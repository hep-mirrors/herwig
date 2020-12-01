  // -*- C++ -*-
  //
  // CMWOneOneOneSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_CMWOneOneOneSplitFn_H
#define HERWIG_CMWOneOneOneSplitFn_H
  //
  // This is the declaration of the CMWOneOneOneSplitFn class.
  //

#include "OneOneOneSplitFn.h"
#include "Herwig/Shower/ShowerAlpha.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"


namespace Herwig {
  
using namespace ThePEG;
/** \ingroup Shower
 *
 * This class provides the concrete implementation 
 * of the CMW enhanced expressions for the
 * splitting function for \f$1\to 11\f$.
 *
 * The kernel uses the same overestimate as the 
 * corresponding OneOneOneSplitFn and thus only needs to 
 * implement the spitting function and ratio to the overestimate.
 *
 * TODO: For a more efficient sampling one needs can rewrite the 
 *       overestimation to contain the alpha_max*Kgmax factors.
 *
 * @see \ref CMWOneOneOneSplitFnInterfaces "The interfaces"
 * defined for CMWOneOneOneSplitFn.
 */
class CMWOneOneOneSplitFn: public Sudakov1to2FormFactor {
    
public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    for(unsigned int ix=0;ix<ids.size();++ix) {
      if(ids[0]->iSpin()!=PDT::Spin1) return false;
    }
    return true;
  }
    
  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   *  Very similar to HalfHalfOneSplitFn.
   *  Here the kernel only contains the soft part multiplied by the
   *  alphas/2pi * Kg from
   *  Nucl.Phys. B349 (1991) 635-654
   *
   */
  virtual double 
  P(const double z, const Energy2 t,
    const IdList & , const bool, const RhoDMatrix &) const {
    Energy2 scale2=t;
    if (!isIS_) scale2 *= pTScale() ? z*(1.-z) : 1.;
    else        scale2 *= pTScale() ? z*(1.-z) :  z;
    return alpha()->value(scale2) * Kg(scale2)/2./Constants::pi/(z*(1.-z));
  }

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
   *  Very similar to HalfHalfOneSplitFn.
   *  Since we use only the 1/1-z part for overestimating the kernel
   *  in the first place we can keep the same overestimation related functions
   *  for the CMW kernels.
   */
  virtual double ratioP(const double z, const Energy2 t,
			const IdList & , const bool, const RhoDMatrix &) const {
    Energy2 scale2=t;
    if (!isIS_) scale2 *= pTScale() ? z*(1.-z) : 1.;
    else        scale2 *= pTScale() ? z*(1.-z) :  z;
    return alpha()->value(scale2)  * Kg(scale2)/2./Constants::pi;
  }

  /**
   * Return the correction term from:
   * Nucl.Phys. B349 (1991) 635-654
   */
  double Kg(Energy2 )const{
    //TODO: Might be t-dependent
    int Nf=5;//alpha_->Nf(t)
    return (3.*(67./18.-1./6.*sqr(Constants::pi))-5./9.*Nf);
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
    assert(PDFfactor==0);
    assert(z>0.&&z<1.);
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
  generatePhiForward(const double z, const Energy2, const IdList &,
		     const RhoDMatrix & rho) {
    assert(rho.iSpin()==PDT::Spin1);
    double modRho = abs(rho(0,2));
    double max = 2.*z*modRho*(1.-z)+sqr(1.-(1.-z)*z)/(z*(1.-z));
    vector<pair<int, Complex> > output;
    output.reserve(3);
    output.push_back(make_pair( 0,(rho(0,0)+rho(2,2))*sqr(1.-(1.-z)*z)/(z*(1.-z))/max));
    output.push_back(make_pair(-2,-rho(0,2)*z*(1.-z)/max));
    output.push_back(make_pair( 2,-rho(2,0)*z*(1.-z)/max));
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
  DecayMEPtr matrixElement(const double z, const Energy2, 
			   const IdList &, const double phi, bool) {
    // calculate the kernal
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
    double omz = 1.-z;
    double root = sqrt(z*omz);
    Complex phase = exp(Complex(0.,1.)*phi);
    (*kernal)(0,0,0) =  phase/root;
    (*kernal)(2,2,2) = -conj((*kernal)(0,0,0));
    (*kernal)(0,0,2) = -sqr(z)/root/phase;
    (*kernal)(2,2,0) = -conj((*kernal)(0,0,2));
    (*kernal)(0,2,0) = -sqr(omz)/root/phase;
    (*kernal)(2,0,2) = -conj((*kernal)(0,2,0));
    (*kernal)(0,2,2) = 0.;
    (*kernal)(2,0,0) = 0.;
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
    
private:

  // Provide information if the kernel is used for initial state.
  // as the pt definition contains an additional factor of z.
  bool isIS_=false;

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

private:
    
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CMWOneOneOneSplitFn & operator=(const CMWOneOneOneSplitFn &) = delete;
  
};
  
}

#endif /* HERWIG_CMWOneOneOneSplitFn_H */
