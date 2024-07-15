// -*- C++ -*-
#ifndef Herwig_ZeroZeroOneEWSplitFn_H
#define Herwig_ZeroZeroOneEWSplitFn_H
//
// This is the declaration of the ZeroZeroOneEWSplitFn class.
//

#include "Sudakov1to2FormFactor.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ZeroZeroOneEWSplitFn class implements the splitting function for
 * \f$\1\to q\1 0\f$ where the spin-1 particles are the W / Z massive
 * electroweak gauge bosons and the spin-0 particle is the massive Higgs
 * boson.
 *
 * @see \ref ZeroZeroOneEWSplitFnInterfaces "The interfaces"
 * defined for ZeroZeroOneEWSplitFn.
 */
class ZeroZeroOneEWSplitFn: public Sudakov1to2FormFactor {

public:

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const {
    if(ids.size()!=3)
      return false;
    if(_couplingValueIm==0.&&_couplingValueRe==0.)
      return false;
    if(ids[0]->iCharge()!=ids[1]->iCharge()+ids[2]->iCharge())
      return false;
    if(ids[0]->iSpin()==PDT::Spin0 && ids[1]->iSpin()==PDT::Spin1 && ids[2]->iSpin()==PDT::Spin0)
      return true;
    else if(ids[0]->iSpin()==PDT::Spin0 && ids[1]->iSpin()==PDT::Spin0 && ids[2]->iSpin()==PDT::Spin1)
      return true;
    else
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
    Complex ghhv(0.,0.);
    getCouplings(ghhv,ids);
    return norm(ghhv)*2*z/(1.-z);
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
			const bool mass, const RhoDMatrix & ) const {
    // ratio in the massless limit
    double val = 1.;
    // the massive limit
    if(mass){
      // get the running mass
      double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
      double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
      double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
      val += m0t2-(1./z)*m1t2+((1.-z)/(4.*z))*m2t2;
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
  virtual double integOverP(const double z, const IdList & ids,
			    unsigned int PDFfactor=0) const {
    assert(PDFfactor==0);
    Complex ghhv(0.,0.);
    getCouplings(ghhv,ids);
    double pre = norm(ghhv)*2.;
    return -pre*(z+log(1.-z));
  }

  /**
   * The concrete implementation of the inverse of the indefinite integral.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double invIntegOverP(const double r, const IdList & ids,
			       unsigned int PDFfactor=0) const {
    assert(PDFfactor==0);
    Complex ghhv(0.,0.);
    getCouplings(ghhv,ids);
    double pre = norm(ghhv)*2.;
    return 1.-exp(-(1.+r/pre));
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
  virtual vector<pair<int,Complex> >
  generatePhiForward(const double , const Energy2 , const IdList & ,
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
  generatePhiBackward(const double , const Energy2 , const IdList & ,
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
				   const IdList & ids, const double phi, bool ) {
    Complex ghhv(0.,0.);
    getCouplings(ghhv,ids);
    // calculate the kernal
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
    double m0t2 = sqr(getParticleData(ids[0]->id())->mass())/t;
    double m1t2 = sqr(getParticleData(ids[1]->id())->mass())/t;
    double m2t2 = sqr(getParticleData(ids[2]->id())->mass())/t;
    Complex phase  = exp(Complex(0.,1.)*phi);
    Complex cphase = conj(phase);
    double sqrtmass = sqrt(m0t2*z*(1.-z)-m1t2*(1.-z)-m2t2*z+z*(1.-z));
    // assign kernel
    (*kernal)(0,0,0) = -phase*ghhv*sqrtmass/(1.-z);        // 111
    (*kernal)(1,0,0) = -sqrt(m2t2)*(1.+z)/sqrt(2.*(1.-z)); // 211 -> 411
    (*kernal)(2,0,0) = cphase*ghhv*sqrtmass/(1.-z);        // 311
    // return the answer
    return kernal;
  }

protected:

  /**
   *   Get the couplings
   */
  void getCouplings(Complex & g, const IdList & ) const  {
    g = Complex(_couplingValueRe,_couplingValueIm);
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
  ZeroZeroOneEWSplitFn & operator=(const ZeroZeroOneEWSplitFn &) = delete;

private:

  /**
   *   numerical value of the splitting coupling to be imported for BSM splittings
   */
  double _couplingValueIm = 0.;
  double _couplingValueRe = 0.;
};

}

#endif /* Herwig_ZeroZeroOneEWSplitFn_H */
