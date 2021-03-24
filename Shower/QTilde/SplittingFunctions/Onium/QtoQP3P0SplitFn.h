// -*- C++ -*-
#ifndef Herwig_QtoQP3P0SplitFn_H
#define Herwig_QtoQP3P0SplitFn_H
//
// This is the declaration of the QtoQP3P0SplitFn class.
//

#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The QtoQP3P0SplitFn class implements the splitting function for \f$q\to q' M_q\bar{q}(^3P_0)'\f$
 *
 * @see \ref QtoQP3P0SplitFnInterfaces "The interfaces"
 * defined for QtoQP3P0SplitFn.
 */
class QtoQP3P0SplitFn: public Sudakov1to2FormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QtoQP3P0SplitFn() : O1_(0.794*GeV*GeV2*GeV2), n_(1), fixedAlphaS_(-1.)
  {}
  //@}

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    // construct the meson PDG code from quark ids and check it
    long id1=ids[0]->id();
    long id2=ids[1]->id();
    long idtest = id1>id2 ? id1*100+id2*10+1 : id2*100+id1*10+1;
    idtest += 10000 + (n_-1)*100000;
    if(abs(ids[2]->id()) != idtest) return false;
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
    return pOver_/z/(1.-z);
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
		const IdList & ids, const bool, const RhoDMatrix &) const {
    Energy m1 = ids[0]->mass();
    Energy M  = m1 + ids[1]->mass();
    double a1 = m1/M;
    double r = sqr(M)/t;
    double W0 = z*sqr( 3.*(-2.+z) +
		       a1*(15.-12.*z+sqr(z) +
			   a1*((-13.+18.*z-5.*sqr(z)) +
			       4.*a1*sqr(1.-z))))/(48.*sqr(a1*sqr(1.-a1*(1.-z))));
    double W1 = (45.-27*z + a1*(6.*(6.-13.*z+5.*sqr(z)) +
				a1*(-407.+591.*z-233.*sqr(z)+9.*z*sqr(z) +
				    a1*(590.-998.*z+498.*sqr(z)-58.*z*sqr(z) +
					8.*a1*(-41.+84.*z-53.*sqr(z)+10.*z*sqr(z) +
					       4.*a1*sqr(1.-z)*(2.-z))))))/(48.*sqr(a1)*sqr(1 - a1*(1 - z)));
    double W2 = (1.-a1)*(-3.+a1*((-22.+12.*z) +
				 a1*((41.-30.*z+3.*sqr(z)) +
				     2.*a1*(-8.+9.*z-sqr(z)))))/(6.*a1*(1 - a1*(1 - z)));
    double W3 = 4./3.*sqr(1.-a1)*a1;
    double ratio = (W0+r*(W1+r*(W2+r*W3)))/pOver_;
    if(ratio>1.) cerr << "ratio greater than 1 in QtoQP3P0SplitFn " << ratio << "\n";
    if(ratio<0.) cerr << "ratio negative       in QtoQP3P0SplitFn " << ratio << " "
		      << r << " " << a1 << " " << z << " " << W0 << " " << W1 << " " << W2 << " " << W3 << "\n";
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
    return pOver_*log(z/(1.-z)); 
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
    return 1./(1.+exp(-r/pOver_));
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
  generatePhiForward(const double , const Energy2 , const IdList &,
		     const RhoDMatrix & ) {
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
  generatePhiBackward(const double , const Energy2, const IdList &,
		      const RhoDMatrix & ) {
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
			   const IdList & ids, const double phi, bool) {
    Energy m1 = ids[0]->mass();
    Energy M  = m1 + ids[1]->mass();
    double a1 = m1/M, a2=1-a1;
    double r = sqr(M)/t;
    Complex ii(0.,1.);
    Complex phase = exp(ii*phi);
    Energy pT = sqrt(z*(1.-z)*t+sqr(M)*(sqr(a1)*z*(1.-z)-sqr(a2)*(1.-z)-z));
    // calculate the kernal N.B. prefactor 1./4./sqrt(3.)/sqrt(z) removed
    DecayMEPtr kernal(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
    (*kernal)(0,0,0) = z*(6.-3.*z-a1*(15.-a1*(13.-5.*z)*(1.-z)+4.*sqr(a1)*sqr(1.-z)-(12.-z)*z))/(a1*sqr(1.-a1*(1.-z)))
      +r*(3.+(1.-2.*a1)*a1*(7.-4.*a1*(1.-z)-5.*z))/a1 - 8.*(1.-a1)*a1*sqr(r)*(1.-a1*(1.-z));
    (*kernal)(1,1,0) = (*kernal)(0,0,0);
    (*kernal)(0,1,0) = double(pT/M)*r/phase*(-8.*(1.-a1)*a1*r+(3.+a1*(7.-2.*a1*(9.-4.*a1*(1.-z)-5.*z)-z))/(a1*(1.-a1*(1.-z))));
    (*kernal)(1,0,0) = -conj((*kernal)(0,1,0));
    return kernal;
  }

protected:
  
  /**
   *  Implementation of the \f$\alpha_S\f$ veto
   */
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
  QtoQP3P0SplitFn & operator=(const QtoQP3P0SplitFn &) = delete;

private:
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy5 O1_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;

  /**
   *  Overestimate of the splitting function
   **/
  static const double pOver_;

  /**
   *  Fixed value of \f$\alpha_S\f$
   */
  double fixedAlphaS_;

};

}

#endif /* Herwig_QtoQP3P0SplitFn_H */
