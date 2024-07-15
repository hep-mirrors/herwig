// -*- C++ -*-
#ifndef Herwig_ZeroZeroZeroEWSplitFn_H
#define Herwig_ZeroZeroZeroEWSplitFn_H
//
// This is the declaration of the ZeroZeroZeroEWSplitFn class.
//

#include "Sudakov1to2FormFactor.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ZeroZeroZeroEWSplitFn class implements the splitting function for
 * \f$\frac12\to q\frac12 h\f$ where the spin-0 higgs particle is a massive scalar boson.
 *
 * @see \ref ZeroZeroZeroEWSplitFnInterfaces "The interfaces"
 * defined for ZeroZeroZeroEWSplitFn.
 */
class ZeroZeroZeroEWSplitFn: public Sudakov1to2FormFactor {

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
    Complex ghhh(0.,0.);
    getCouplings(ghhh,ids);
    return norm(ghhh)/(2.*z*(1.-z));
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
  virtual double ratioP(const double z, const Energy2 t, const IdList & ,
			const bool , const RhoDMatrix & ) const {
    return z*(1.-z)/t*GeV2;
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
  generatePhiForward(const double, const Energy2, const IdList &,
	      const RhoDMatrix &) {
    // scalar so no dependence
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
  virtual vector<pair<int,Complex> >
  generatePhiBackward(const double, const Energy2, const IdList &,
		      const RhoDMatrix &) {
    // scalar so no dependence
    assert(false);
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
   *   Get the couplings
   * @param g The couplings ( for the SM case, this is g_SM/(e*mH^2) )
   * @param ids The PDG codes for the particles in the splitting.
   */
  void getCouplings(Complex & g, const IdList & ids) const;

  /**
   *   Get the couplings with running masses
   * @param g The couplings ( for the SM case, this is g_SM/(e*mH^2) )
   * @param ids The PDG codes for the particles in the splitting.
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   */
  void getCouplings(Complex & g, const IdList & ids, const Energy2 t) const;

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
  ZeroZeroZeroEWSplitFn & operator=(const ZeroZeroZeroEWSplitFn &) = delete;

private:

  /**
   *  1 / sin theta_w
   */
  double gw_;

  /**
   *   numerical value of the splitting coupling to be imported for BSM splittings
   */
  double _couplingValueIm = 0.;
  double _couplingValueRe = 0.;

  /**
   * Pointer to the SM object.
   */
  tcHwSMPtr _theSM;

};

}

#endif /* Herwig_ZeroZeroZeroEWSplitFn_H */
