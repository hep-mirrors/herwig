// -*- C++ -*-
//
// OneHalfHalfSplitFn.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_OneHalfHalfSplitFn_H
#define HERWIG_OneHalfHalfSplitFn_H
//
// This is the declaration of the OneHalfHalfSplitFn class.
//

#include "Sudakov1to2FormFactor.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This class provides the concrete implementation of the exact leading-order
 * splitting function for \f$1\to \frac12\frac12\f$. 
 *
 * In this case the splitting function is given by
 * \f[P(z,t) =C\left(1-2z(1-z)+2\frac{m_q^2}{t}\right),\f]
 * where \f$C\f$ is the corresponding colour factor
 * Our choice for the overestimate is 
 * \f[P_{\rm over}(z) = C,\f]
 * therefore the integral is
 * \f[\int P_{\rm over}(z) {\rm d}z =Cz,\f]
 * and its inverse is
 * \f[\frac{r}{C}\f]
 *
 * @see \ref OneHalfHalfSplitFnInterfaces "The interfaces"
 * defined for OneHalfHalfSplitFn.
 */
class OneHalfHalfSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    if(ids[1]!=ids[2]->CC()) return false;
    if(ids[1]->iSpin()!=PDT::Spin1Half) return false;
    if(ids[0]->iSpin()!=PDT::Spin1) return false;
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
  virtual double overestimateP(const double , const IdList & ) const {
    return 1.;
  }

  /**
   * The concrete implementation of the
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   * @param rho The spin density matrix
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			const bool mass, const RhoDMatrix &) const {
    double zz = z*(1.-z);
    double val = 1.-2.*zz;
    if(mass) {
      Energy m = ids[1]->mass();
      val+= 2.*sqr(m)/t;
    }
    return val;
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
  virtual double integOverP(const double z,  const IdList & ids, 
			    unsigned int PDFfactor=0) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r,  const IdList & ids, 
			       unsigned int PDFfactor=0) const;
  //@}

  /**
   * Method to calculate the azimuthal angle
   * @param particle The particle which is branching
   * @param showerkin The ShowerKinematics object
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
   * Method to calculate the azimuthal angle
   * @param particle The particle which is branching
   * @param showerkin The ShowerKinematics object
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> >
  generatePhiBackward(const double, const Energy2, const IdList &,
		      const RhoDMatrix &) { 
    // no dependance
    return {{ {0, 1.} }};
  }
  
  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   */
  virtual DecayMEPtr matrixElement(const double z, const Energy2 t, 
                                   const IdList & ids, const double phi, bool timeLike);

protected:
  
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
   *   Interpolation routine
   */
  void findZero();

  /**
   *   Override PDF veto function
   */
  virtual double PDFVetoRatio(const Energy2 t, const double x, const double z,
                              const tcPDPtr parton0, const tcPDPtr parton1, double factor) const;

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
  OneHalfHalfSplitFn & operator=(const OneHalfHalfSplitFn &) = delete;

private:

  /**
   *   Storage of masses for sampling of q2
   */
  map<cPDFPtr,vector<Energy2>> mq_;

  /**
   *  Current mass
   */
  Energy2 mq2_;
};

}

#endif /* HERWIG_OneHalfHalfSplitFn_H */
