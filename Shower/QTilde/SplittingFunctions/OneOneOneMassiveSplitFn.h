// -*- C++ -*-
//
// OneOneOneSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_OneOneOneMassiveSplitFn_H
#define HERWIG_OneOneOneMassiveSplitFn_H
//
// This is the declaration of the OneOneOneMassiveSplitFn class.
//

#include "Sudakov1to2FormFactor.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$1\to 11\f$ where the emitting particle is massi e
 *
 * In this case the splitting function is given by
 * \f[P(z) =2C\left( \frac{z}{1-z}-\frac{m^2}{q^2-m^2} +2\rho_{00}\frac{(1-z)^2m^2}{q^2-m^2} + (1-\rho_{00})\left(\frac{1-z}{z}+z(1-z)-\frac{(1-z)^2m^2}{q^2-m^2}\right)\right),\f]
 * where \f$C=\f$ is the corresponding colour factor.
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = 2C\left(\frac1z+\frac1{1-z}\right),\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z =2C\ln\left(\frac{z}{1-z}\right),\f]
 * and its inverse is
 * \f[\frac{\exp\left(\frac{r}{2C}\right)}{(1+\exp\left(\frac{r}{2C}\right)}\f]
 *
 *
 * @see \ref OneOneOneMassiveSplitFnInterfaces "The interfaces"
 * defined for OneOneOneMassiveSplitFn.
 */
class OneOneOneMassiveSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    for(unsigned int ix=0;ix<ids.size();++ix) {
      if(ids[0]->iSpin()!=PDT::Spin1) return false;
    }
    return true;
  }

  /**
   * The concrete implementation of the overestimate of the splitting function,
   * \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList &) const {
    return 2.*(1/z + 1/(1.-z)); 
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
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			const bool, const RhoDMatrix & rho) const {
    Energy2 m2 = sqr(ids[0]->mass());
    double rho00 = rho(1,1).real();
    return (sqr(z)-z*(1.-z)*m2/t
            + (1.-rho00)*sqr(1.-z)*( 1. + sqr(z) -z*(1.-z)*m2/t )
            + 2.*rho00*z*(1.-z)*sqr(1.-z)*m2/t);
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
  virtual double integOverP(const double z, const IdList &,
			    unsigned int PDFfactor=0) const {
    assert(PDFfactor==0);
    assert(z>0.&&z<1.);
    return 2.*log(z/(1.-z));
  }

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r, const IdList &,
			       unsigned int PDFfactor=0) const {
    assert(PDFfactor==0);
    return 1./(1.+exp(-0.5*r)); 
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
  virtual vector<pair<int,Complex> >
  generatePhiForward(const double z, const Energy2 t, const IdList & ids,
		     const RhoDMatrix &);

  /**
   * Method to calculate the azimuthal angle for backward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiBackward(const double z, const Energy2 t, const IdList & ids,
		      const RhoDMatrix &);
  
  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   */
  virtual DecayMEPtr matrixElement(const double z, const Energy2 t, 
				   const IdList & ids, const double phi, bool timeLike);

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
  OneOneOneMassiveSplitFn & operator=(const OneOneOneMassiveSplitFn &) = delete;

};

}

#endif /* HERWIG_OneOneOneMassiveSplitFn_H */
