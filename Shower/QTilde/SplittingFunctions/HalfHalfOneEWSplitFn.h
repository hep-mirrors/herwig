// -*- C++ -*-
#ifndef Herwig_HalfHalfOneEWSplitFn_H
#define Herwig_HalfHalfOneEWSplitFn_H
//
// This is the declaration of the HalfHalfOneEWSplitFn class.
//

#include "Sudakov1to2FormFactor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The HalfHalfOneEWSplitFn class implements the splitting function for 
 * \f$\frac12\to q\frac12 1\f$ where the spin-1 particle is a massive electroweak gauge boson.
 *
 * @see \ref HalfHalfOneEWSplitFnInterfaces "The interfaces"
 * defined for HalfHalfOneEWSplitFn.
 */
class HalfHalfOneEWSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const;

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
  virtual double overestimateP(const double z, const IdList & ids) const {
    Complex gL(0.,0.),gR(0.,0.);
    getCouplings(gL,gR,ids);
    return 2.*max(norm(gL),norm(gR))/(1.-z); //FIXME//
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
			const bool mass, const RhoDMatrix & rho) const;

  /**
   * The concrete implementation of the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double integOverP(const double z, const IdList & ids, 
			    unsigned int PDFfactor=0) const;

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */ 
  virtual double invIntegOverP(const double r, const IdList & ids, 
			       unsigned int PDFfactor=0) const;
  //@}

  /**
   * Method to calculate the azimuthal angle
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> >
  generatePhiForward(const double, const Energy2, const IdList &, const RhoDMatrix &) {
    // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
    // and rest = splitting function for Tr(rho)=1 as required by defn
    return vector<pair<int, Complex> >(1,make_pair(0,1.));
  }

  /**
   * Method to calculate the azimuthal angle for backward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiBackward(const double, const Energy2, const IdList &, const RhoDMatrix &) {
    // no dependence on the spin density matrix, dependence on off-diagonal terms cancels
    // and rest = splitting function for Tr(rho)=1 as required by defn
    return vector<pair<int, Complex> >(1,make_pair(0,1.));
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
   *   Get the couplings
   */
  void getCouplings(Complex & gL, Complex & gR, const IdList & ids) const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HalfHalfOneEWSplitFn & operator=(const HalfHalfOneEWSplitFn &) = delete;

private:

  /**
   *  Z couplings
   */
  map<long,pair<double,double> > gZ_;

  /**
   *  W couplings
   */
  double gWL_;

  /** 
   *   numerical value of the splitting coupling to be imported for BSM splittings
   */
  double _couplingValueLeftRe = 0.; 
  double _couplingValueLeftIm = 0.; 
  double _couplingValueRightRe = 0.; 
  double _couplingValueRightIm = 0.;

};

}

#endif /* Herwig_HalfHalfOneEWSplitFn_H */
