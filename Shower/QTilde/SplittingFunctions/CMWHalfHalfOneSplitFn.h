  // -*- C++ -*-
  //
  // CMWHalfHalfOneSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef HERWIG_CMWHalfHalfOneSplitFn_H
#define HERWIG_CMWHalfHalfOneSplitFn_H
  //
  // This is the declaration of the CMWHalfHalfOneSplitFn class.
  //

#include "HalfHalfOneSplitFn.h"
#include "Herwig/Shower/ShowerAlpha.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"


namespace Herwig {
  
using namespace ThePEG;
  
/** \ingroup Shower
 *
 * This class provides the concrete implementation
 * of the CMW enhanced expressions for the
 * splitting function for \f$\frac12\to q\frac12 1\f$.
 *
 * The kernel uses the same overestimate as the
 * corresponding HalfHalfOneSplitFn and thus only needs to
 * implement the spitting function and ratio to the overestimate.
 *
 * TODO: For a more efficient sampling one needs can rewrite the
 *       overestimation to contain the alpha()max*Kgmax factors.
 *
 * @see \ref CMWHalfHalfOneSplitFnInterfaces "The interfaces"
 * defined for CMWHalfHalfOneSplitFn.
 */
class CMWHalfHalfOneSplitFn: public Sudakov1to2FormFactor {
    
public:
    
  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  bool accept(const IdList & ids) const {
    // 3 particles and in and out fermion same
    if(ids.size()!=3 || ids[0]!=ids[1]) return false;
    if(ids[0]->iSpin()!=PDT::Spin1Half ||
       ids[2]->iSpin()!=PDT::Spin1) return false;
    return true;
  }

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  double overestimateP(const double z, const IdList & ) const { 
    return 2./(1.-z); 
  }

  /**
   *  Very similar to HalfHalfOneSplitFn.
   *  Since we use only the 1/1-z part for overestimating the kernel 
   *  in the first place we can keep the same overestimation related functions
   *  for the CMW kernels.
   */
  virtual double ratioP(const double z, const Energy2 t,
                                     const IdList &, const bool , const RhoDMatrix & ) const {
    Energy2 scale2=t;
    //See pt definitions in QTildeSudakov.cc
    // Note: t here is t * f(z)
    if (!isIS_)    scale2 *= pTScale() ? z*(1.-z): 1.;
    else           scale2 *= pTScale() ? z*(1.-z): z ;
    return 0.5 * Kg(scale2) * alpha()->value(scale2)/2./Constants::pi;
  }
    
  /**
   * Return the correction term from:
   * Nucl.Phys. B349 (1991) 635-654
   */
  double Kg(Energy2 ) const{
    //TODO: Should be t-dependent
    int Nf=5;//alpha()->Nf(t)
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
  double integOverP(const double z, const IdList &, 
		    unsigned int PDFfactor=0) const {
    switch (PDFfactor) {
    case 0:
      return -2.*Math::log1m(z);
    case 1:
      return  2.*log(z/(1.-z));
    case 2:
      return  2./(1.-z);
    case 3:
    default:
      throw Exception() << "CMWHalfHalfOneSplitFn::integOverP() invalid PDFfactor = "
			<< PDFfactor << Exception::runerror;
    }
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
		       unsigned int PDFfactor) const {
    switch (PDFfactor) {
    case 0:
      return 1. - exp(- 0.5*r); 
    case 1:
      return 1./(1.-exp(-0.5*r));
    case 2:
      return 1.-2./r;
    case 3:
    default:
      throw Exception() << "CMWHalfHalfOneSplitFn::invIntegOverP() invalid PDFfactor = "
			<< PDFfactor << Exception::runerror;
    } 
  }
  //@}

  /**
   * Method to calculate the azimuthal angle
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  vector<pair<int,Complex> >
  generatePhiForward(const double, const Energy2, const IdList & , const RhoDMatrix &) {
    // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
    // and rest = splitting function for Tr(rho)=1 as required by defn
    return {{ {0, 1.} }};
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
  generatePhiBackward(const double, const Energy2, const IdList & , const RhoDMatrix &) {
    // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
    // and rest = splitting function for Tr(rho)=1 as required by defn
    return {{ {0, 1.} }};
  }
  
  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   */
  DecayMEPtr matrixElement(const double z, const Energy2 t, 
						const IdList & ids, const double phi,
						bool timeLike) {
    // calculate the kernal
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
    Energy m = !timeLike ? ZERO : ids[0]->mass();
    double mt = m/sqrt(t);
    double root = sqrt(1.-(1.-z)*sqr(m)/z/t);
    double romz = sqrt(1.-z); 
    double rz   = sqrt(z);
    Complex phase = exp(Complex(0.,1.)*phi);
    (*kernal)(0,0,0) = -root/romz*phase;
    (*kernal)(1,1,2) =  -conj((*kernal)(0,0,0));
    (*kernal)(0,0,2) =  root/romz*z/phase;
    (*kernal)(1,1,0) = -conj((*kernal)(0,0,2));
    (*kernal)(1,0,2) =  mt*(1.-z)/rz;
    (*kernal)(0,1,0) =  conj((*kernal)(1,0,2));
    (*kernal)(0,1,2) =  0.;
    (*kernal)(1,0,0) =  0.;
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
  CMWHalfHalfOneSplitFn & operator=(const CMWHalfHalfOneSplitFn &) = delete;
  
};
}

#endif /* HERWIG_CMWHalfHalfOneSplitFn_H */
