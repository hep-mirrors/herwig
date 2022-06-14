// -*- C++ -*-
#ifndef Herwig_KMatrix_H
#define Herwig_KMatrix_H
//
// This is the declaration of the KMatrix class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Decay/IsoSpin.h"
#include "KMatrix.fh"
#include <boost/numeric/ublas/matrix.hpp>

namespace Herwig {
namespace ublas = boost::numeric::ublas;
using namespace ThePEG;

/**
 * The KMatrix class is a base class for the implementation of
 * K-matrix parameterizations in Herwig
 *
 * @see \ref KMatrixInterfaces "The interfaces"
 * defined for KMatrix.
 */
class KMatrix: public Interfaced {

public:
  
  /**
   * Enum for the possible channels
   */
  enum Channels { PiPi, KPi, KEta, KEtaPrime, KK, EtaEta,EtaEtaPrime,FourPi};

public:

  /**
   * The default constructor.
   */
  KMatrix(FlavourInfo flavour=FlavourInfo(),
	  vector<Channels> channels=vector<Channels>(),
	  vector<Energy2> poles=vector<Energy2>(),
	  vector<vector<Energy> > g=vector<vector<Energy> >());

  /**
   *   The quantum numbers of the K-matrix
   */
  FlavourInfo flavourInfo() const {
    return flavour_;
  };

  /**
   * Compute the K-matrix for a given scale
   * @param s The scale
   * @param Whether or not to multiply by \f$\prod_i(1-s/m^2_i)\f$ to regularise the poles
   */
  virtual ublas::matrix<double> K(Energy2 s, bool multiplyByPoles=false) const = 0;

  /**
   *   The \f$\rho\f$ matrix
   */
  virtual ublas::matrix<Complex> rho(Energy2 s) const;
  
  /**
   *  Vector containing the locations of the poles
   */
  const vector<Energy2> & poles() const {return poles_;}

  /**
   *  Access the couplings of the poles
   */
  const vector<vector<Energy> > & poleCouplings() const {return g_;}

  /**
   * Compute the amplitdes given the \f$P\f$-vector
   * @param Whether or not to multiply by \f$\prod_i(1-s/m^2_i)\f$ to regularise the poles
   */
  virtual ublas::vector<Complex>
  amplitudes(Energy2 s, ublas::vector<Complex> pVector,
	     bool multiplyByPoles=false) const;

  /**
   *  The number of channels
   */
  unsigned int numberOfChannels() {return channels_.size();}
  
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
  KMatrix & operator=(const KMatrix &) = delete;

private:

  /**
   *   The quantum numbers for the K-matrix
   */
  FlavourInfo flavour_;
  
  /**
   *   The mesons in the various channels
   */
  vector<Channels> channels_;

  /**
   *  The positions of the poles
   */
  vector<Energy2> poles_;

  /**
   *  Couplings for the resonances
   */
  vector<vector<Energy> > g_;

private:

  /**
   *  Common masses for the \f$\rho\f$ matrix
   */
  //@{
  /**
   *   The charged pion mass
   */
  Energy mPiPlus_;

  /**
   *   The neutral pion mass
   */
  Energy mPi0_;

  /**
   *   The charged kaon mass
   */
  Energy mKPlus_;

  /**
   *   The neutral kaon mass
   */
  Energy mK0_;

  /**
   *  The \f$\eta\f$ mass
   */
  Energy mEta_;

  /**
   *  The \f$\eta^\prime\f$ mass
   */
  Energy mEtaPrime_;
  //@}
};

}

#endif /* Herwig_KMatrix_H */
