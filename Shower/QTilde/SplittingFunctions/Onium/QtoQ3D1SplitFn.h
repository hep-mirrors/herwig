// -*- C++ -*-
#ifndef Herwig_QtoQ3D1SplitFn_H
#define Herwig_QtoQ3D1SplitFn_H
//
// This is the declaration of the QtoQ3D1SplitFn class.
//

#include "Herwig/Shower/QTilde/SplittingFunctions/Sudakov1to2FormFactor.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The QtoQ3D1SplitFn class implements the splitting function for \f$q\to q M_q\bar{q}(^3D_1)\f$.
 *
 * @see \ref QtoQ3D1SplitFnInterfaces "The interfaces"
 * defined for QtoQ3D1SplitFn.
 */
class QtoQ3D1SplitFn: public Sudakov1to2FormFactor {

public:

  /**
   * The default constructor.
   */
  QtoQ3D1SplitFn() : O1_(0.131*GeV*pow<3,1>(GeV2)), n_(1), fixedAlphaS_(-1.)
  {}

  /**
   *  Concrete implementation of the method to determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  bool accept(const IdList & ids) const {
    if(ids.size()!=3) return false;
    // construct the meson PDG code from quark ids and check it
    long idtest = ids[0]->id()*110+3 + 30000 + (n_-1)*100000;
    if(ids[2]->id() != idtest) return false;
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
		const IdList & ids, const bool, const RhoDMatrix &) const;
  
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
			   const IdList & ids, const double phi, bool);

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
  QtoQ3D1SplitFn & operator=(const QtoQ3D1SplitFn &) = delete;

private:
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy7 O1_;

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

#endif /* Herwig_QtoQ3D1SplitFn_H */
