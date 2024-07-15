// -*- C++ -*-
#ifndef Herwig_DalitzKMatrix_H
#define Herwig_DalitzKMatrix_H
//
// This is the declaration of the DalitzKMatrix class.
//

#include "DalitzResonance.h"
#include "Herwig/Decay/FormFactors/KMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DalitzKMatrix class allows the use of \f$K\f$-matrices in Dalitz decays
 */
class DalitzKMatrix: public DalitzResonance {

public:

  /**
   * The default constructor.
   */
  DalitzKMatrix()
  {}

  /**
   *  Constructor specifiying the parameters
   */
  DalitzKMatrix(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		unsigned int d1, unsigned int d2, unsigned int s,
		double mag, double phi, InvEnergy rr,
		unsigned int imat, unsigned int chan,
		Energy2 sc, unsigned int itype,
		vector<pair<double,double> > beta,
		vector<pair<double,vector<double > > > coeffs)
    : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr),
      imat_(imat), channel_(chan), sc_(sc), expType_(itype), coeffs_(coeffs) {
    beta_.clear();
    for(unsigned int ix=0;ix<beta.size();++ix) {
      beta_.push_back(beta[ix].first*exp(Complex(0.,beta[ix].second)));
    }
  }

public:

  /**
   *  Return the Breit-Wigner times the form factor
   */
  virtual Complex BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const;

  /**
   *  Output the parameters
   */
  virtual void dataBaseOutput(ofstream & output);

public:

  /**
   *  Set the K-matrix
   */
  void setKMatrix(KMatrixPtr mat) {kMatrix_ = mat;}

  /**
   *  Location of the matrix
   */
  unsigned int imatrix() const {return imat_;}

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

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DalitzKMatrix & operator=(const DalitzKMatrix &) = delete;

private:

  /**
   *   The K-matrix for the channel
   */
  KMatrixPtr kMatrix_;

  /**
   *  Which \f$K-matrix\f$ to do
   */
  unsigned int imat_ = 0;

  /**
   *  Which channel to use from the K-matrix
   */
  unsigned int channel_ = 0;

  /**
   *  Expansion point for the constant terms
   */
  Energy2 sc_ = ZERO;

  /**
   *  Coefficients of the poles
   */
  vector<Complex> beta_;

  /**
   *  Type of expansion
   */
  unsigned int expType_ = 0;

  /**
   *  Coefficients for the series expansion
   */
  vector<pair<double,vector<double > > > coeffs_;

};

}

#endif /* Herwig_DalitzKMatrix_H */
