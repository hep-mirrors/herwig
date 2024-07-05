// -*- C++ -*-
#ifndef Herwig_HalfHalfZeroEWSplitFn_H
#define Herwig_HalfHalfZeroEWSplitFn_H
//
// This is the declaration of the HalfHalfZeroEWSplitFn class.
//

#include "Sudakov1to2FormFactor.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The HalfHalfZeroEWSplitFn class implements the splitting function for
 * \f$\frac12\to q\frac12 h\f$ where the spin-0 higgs particle is a massive scalar boson.
 *
 * @see \ref HalfHalfZeroEWSplitFnInterfaces "The interfaces"
 * defined for HalfHalfZeroEWSplitFn.
 */
class HalfHalfZeroEWSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    if(ids[2]->iSpin()==PDT::Spin0 && _couplingValue0Re==0 && _couplingValue0Im==0 && _couplingValue1Re==0 && _couplingValue1Im==0) {
      if(ids[0]->id()==ids[1]->id() && (ids[0]->id()==4 || ids[0]->id()==5 || ids[0]->id()==6)) return true;
    }
    else if(ids[2]->iSpin()==PDT::Spin0 && !(_couplingValue0Re==0 && _couplingValue0Im==0 && _couplingValue1Re==0 && _couplingValue1Im==0)) {
      if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge()) return false;
      if((abs(ids[0]->id())>=1 && abs(ids[0]->id())<=6) && (abs(ids[1]->id())>=1 && abs(ids[1]->id())<=6)) return true;
    }
    return false;
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
  virtual double overestimateP(const double z, const IdList & ids) const {
    Complex gH0(0.,0.);
    Complex gH1(0.,0.);
    getCouplings(gH0,gH1,ids);
    return (norm(gH0)+norm(gH1)+2*abs((gH0*gH1).real())+2*abs((gH0*gH1).imag()))*(1.-z)/2.;
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
  generatePhiForward(const double, const Energy2, const IdList &,
                     const RhoDMatrix &) {
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
  generatePhiBackward(const double, const Energy2, const IdList & ,
		      const RhoDMatrix &) {
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
   *   Get the couplings without running masses
   */
  void getCouplings(Complex & gH0, Complex & gH1, const IdList & ids) const;

  /**
   *   Get the couplings with running masses
   */
  void getCouplings(Complex & gH0, Complex & gH1, const IdList & ids, const Energy2 t) const;

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
  HalfHalfZeroEWSplitFn & operator=(const HalfHalfZeroEWSplitFn &) = delete;

private:

  /**
   *  Higgs couplings
   */
  double ghqq_;

  /**
   *   numerical value of the splitting coupling to be imported for BSM splittings
   */
  double _couplingValue0Re = 0.;
  double _couplingValue0Im = 0.;
  double _couplingValue1Re = 0.;
  double _couplingValue1Im = 0.;

  /**
   * Pointer to the SM object.
   */
  tcHwSMPtr _theSM;

};

}

#endif /* Herwig_HalfHalfZeroEWSplitFn_H */
