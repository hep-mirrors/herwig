// -*- C++ -*-
//
// HalfHalfOneSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HalfHalfOneSplitFn_H
#define HERWIG_HalfHalfOneSplitFn_H
//
// This is the declaration of the HalfHalfOneSplitFn class.
//

#include "Sudakov1to2FormFactor.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$\frac12\to q\frac12 1\f$. 
 *
 *  In this case the splitting function is given by
 * \f[P(z,t) =C\left(\frac{1+z^2}{1-z}-2\frac{m^2_q}{t}\right),\f]
 * where \f$C\f$ is the corresponding colour factor.
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = \frac{2C}{1-z},\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z = -2C\ln(1-z),\f]
 * and its inverse is
 * \f[1-\exp\left(\frac{r}{2C}\right).\f]
 *
 * @see \ref HalfHalfOneSplitFnInterfaces "The interfaces"
 * defined for HalfHalfOneSplitFn.
 */
class HalfHalfOneSplitFn: public Sudakov1to2FormFactor {

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
                const IdList & ids, const bool mass, const RhoDMatrix & ) const {
    double val = 1. + sqr(z);
    if(mass) {
      Energy m = ids[0]->mass();
      val -= 2.*sqr(m)*(1.-z)/t;
    } 
    return 0.5*val;
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
                    unsigned int PDFfactor) const {
    switch (PDFfactor) {
    case 0:
      return -2.*Math::log1m(z);
    case 1:
      return  2.*log(z/(1.-z));
    case 2:
      return  2./(1.-z);
    case 3:
    default:
      throw Exception() << "HalfHalfOneSplitFn::integOverP() invalid PDFfactor = "
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
      throw Exception() << "HalfHalfOneSplitFn::invIntegOverP() invalid PDFfactor = "
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

public:

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HalfHalfOneSplitFn & operator=(const HalfHalfOneSplitFn &) = delete;

};

}

#endif /* HERWIG_HalfHalfOneSplitFn_H */
