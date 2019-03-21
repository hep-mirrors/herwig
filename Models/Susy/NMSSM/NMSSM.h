// -*- C++ -*-
#ifndef HERWIG_NMSSM_H
#define HERWIG_NMSSM_H
//
// This is the declaration of the NMSSM class.
//

#include "Herwig/Models/Susy/MSSM.h"
#include "NMSSM.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * Here is the documentation of the NMSSM class.
 *
 * @see \ref NMSSMInterfaces "The interfaces"
 * defined for NMSSM.
 */
class NMSSM: public MSSM {

public:

  /**
   * The default constructor.
   */
  NMSSM() : _lambda(0.), _kappa(0.), _theAlambda(0.*MeV), 
	    _theAkappa(0.*MeV), _lambdaVEV(0.*MeV),
		_MQ3(0.*MeV), _MU2(0.*MeV) 
  {}

public:

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  double lambda() const {
    return _lambda;
  }

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  double kappa() const {
    return _kappa;
  }

  /**
   *  The V.E.V of the extra singlet field scaled
   * by \f$ lambda\f$, 
   */
  Energy lambdaVEV() const {
    return _lambdaVEV;
  }
  
  /**
   * Soft trilinear \f$SH_2 H_1\f$ coupling
   */
  Energy trilinearLambda() const {
    return _theAlambda;
  }

  /**
   * Soft cubic \f$S\f$ coupling
   */
  Energy trilinearKappa() const {
    return _theAkappa;
  }
      /**
   *  left 3rd generation scalar quark mass
   */
  Energy MQ3() const {
    return _MQ3;
  }

  /**
   * right scalar top mass
   */
  Energy MU2() const {
    return _MU2;
  }
  //@}

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

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

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
  NMSSM & operator=(const NMSSM &) = delete;

private:

  /**
   *  The NMSSM couplings
   */
  //@{
  /**
   *  Superpotential \f$\lambda\f$ term
   */
  double _lambda;

  /**
   *  Superpotential \f$\kappa\f$ coupling
   */
  double _kappa;

  /**
   * Soft trilinear \f$SH_2 H_1\f$ coupling
   */
  Energy _theAlambda;

  /**
   * Soft cubic \f$S\f$ coupling
   */
  Energy _theAkappa;

  /**
   *  The V.E.V of the extra singlet field scaled
   * by \f$ lambda\f$
   */
  Energy _lambdaVEV;
      /**
   * left 3rd generation scalar quark mass
   */
  Energy _MQ3;

  /**
   *  right scalar top mass
   */
  Energy _MU2;
  //@}
};

}

#endif /* HERWIG_NMSSM_H */
