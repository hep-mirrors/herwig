// -*- C++ -*-
//
// UEDBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDBase_H
#define HERWIG_UEDBase_H
//
// This is the declaration of the UEDBase class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"

#include "UEDBase.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class serves as a base class for all UED models. It stores the
 * values of the inverse radius and the product \f$\Lambda R \f$  and has functions
 * to calculate the radiative corrections to the nth level KK excitations.
 *
 * To use this class for n > 1 simply inherit off it, calculate the necessary masses
 * using the provided functions for the new excitations and add the
 * appropriate vertices.
 *
 * @see \ref UEDBaseInterfaces "The interfaces"
 * defined for UEDBase.
 */
class UEDBase: public BSMModel {

public:

  /** Typedef for ID-Mass pair. */
  typedef pair<long, Energy> IDMassPair;

  /** Typedef for unsigned int/double map to store Weinburg angles.*/
  typedef map<unsigned int, double> WAMap;

public:

  /**
   * The default constructor.
   */
  UEDBase();

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

public:

  /** @name Public Access Functions.*/
  //@{
  /**
   * Return the compactification radius
   */
  InvEnergy compactRadius() const {
    return 1./theInvRadius;
  }

  /**
   * Return the Weinburg mixing angle for any level.
   */
  double sinThetaN(const unsigned int n) const;

  /**
   * Return the Weinburg mixing angle for \f$n = 1\f$
   */
  double sinThetaOne() const {
    return theSinThetaOne;
  }
  //@}

protected:

  /**
   * Add a new ID,mass pair to the mass storage
   * @param elem The element to add in to storage
   */
  void addMassElement(IDMassPair elem) {
    theMasses.push_back(elem);
  }

  /**
   * Add a new mixing angle to the storage
   * @param n The level
   * @param val The value
   */
  void addMixingAngle(const unsigned int n, 
		      const double val) {
    theMixingAngles.insert(make_pair(n, val));
  }
  
private:

  /** @name Utility Functions for calculating masses. */
  //@{
  /**
   * Calculate the radiative corrections to the masses of the KK excitations
   * @param n The KK-level for which to calculate the masses. 
   */
  void calculateKKMasses(const unsigned int n);

  /**
   * Calculate the radiative corrections to the spin-0 and spin-1 
   * masses of the KK excitations
   * @param n The KK-level for which to calculate the masses. 
   */
  void bosonMasses(const unsigned int n);
  
  /**
   * Calculate the radiative corrections to the spin-1/2
   * masses of the KK excitations.
   * @param n The KK-level for which to calculate the masses.
   */
  void fermionMasses(const unsigned int n);

  /**
   * Reset the mass of the ParticleData object
   *@param id The id of the particles mass to reset
   *@param value The new mass
   */  
  void resetMass(long id, Energy value);

  /**
   * Calculate the Weinburg Mixing angle for the appropriate level.
   * @param n The KK-level for which to calculate the mixing angle.
   */
  double calculateMixingAngle(const unsigned int n);
  //@}
  
  /**
   * Write out a spectrum file ordered in mass (name can be set by an interface).
   */
  void writeSpectrum();

  /**
   * A predicate for sorting the list of masses.
   */
  static bool lowerMass(const pair<long, Energy> & p1, 
			const pair<long, Energy> & p2) {
    return p1.second < p2.second;
  }
  
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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDBase & operator=(const UEDBase &) = delete;

private:
  
  /**
   * Whether to calculate the radiative corrections to the KK masses
   */
  bool theRadCorr;
  
  /**
   * Store the radius of the compactified dimension.
   */
  Energy theInvRadius;
  
  /**
   * The value of \f$\Lambda R \f$.
   */
  double theLambdaR;

  /**
   * The boundary mass term for the Higgs.
   */
  Energy theMbarH;

  /**
   * The values of \f$\sin\theta_N\f$
   */
  WAMap theMixingAngles;

  /**
   * Store \f$\sin\theta_1\f$ for faster access
   */
  double theSinThetaOne;
  
  /**
   * Store the masses of the new particles
   */
  vector<IDMassPair> theMasses;
  
  /**
   * The value of the vacuum expectation value of the higgs field.
   */
  Energy theVeV;

  /**
   *  Include SM masses in calculation of KK masses
   */
  bool includeSMMass_;

  /**
   *  Use fixed couplings for the mass calculation
   */
  bool fixedCouplings_;

  /**
   *  Include gauge boson mixing
   */
  bool includeGaugeMixing_;
  
  /** @name The level 1 UED vertices. */
  //@{
  /**
   * The \f$\bar{f}^{(1)}f^{(1)}Z^{(0)}\f$
   */
  AbstractFFVVertexPtr theF1F1Z0Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(1)}g^{(0)}\f$
   */
  AbstractFFVVertexPtr theF1F1G0Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(0)}g^{(1)}\f$
   */
  AbstractFFVVertexPtr theF1F0G1Vertex;

  /**
   * The \f$g^{(1)}g^{(1)}g\f$ vertex
   */
  AbstractVVVVertexPtr theG1G1G0Vertex;

  /**
   * The \f$g\,g\,g^{(1)},g^{(1)}\f$ vertex
   */
  AbstractVVVVVertexPtr theG0G0G1G1Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(1)}\gamma\f$
   */
  AbstractFFVVertexPtr theF1F1P0Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(1)}W\f$
   */
  AbstractFFVVertexPtr theF1F1W0Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(0)}W^{(1)}\f$
   */
  AbstractFFVVertexPtr theF1F0W1Vertex;

  /**
   * The \f$\bar{f}^{(1)}f^{(0)}H^{(1)}\f$
   */
  AbstractFFSVertexPtr theF1F0H1Vertex;

  /**
   * The \f$ A^\mu_{(0)}H^+_{(1)}H-_{(1)}\f$
   */
  AbstractVSSVertexPtr theP0H1H1Vertex;

  /**
   * The \f$ Z^\mu_{(0)}H^+_{(1)}H-_{(1)}\f$
   */
  AbstractVSSVertexPtr theZ0H1H1Vertex;

  /**
   * The \f$ W^\pm_{\mu(0)}A_{(1)}H^\mp_{(1)}\f$
   */
  AbstractVSSVertexPtr theW0A1H1Vertex;

  /**
   * The \f$ Z^\mu_{\mu(0)}A_{(1)}h_{(1)}\f$
   */
  AbstractVSSVertexPtr theZ0A1h1Vertex;
  
  /**
   * The \f$W^{(1)}Z^{(1)}W_{(0)}\f$ vertex
   */
  AbstractVVVVertexPtr theW0W1W1Vertex;
  //@}
};


}

#endif /* HERWIG_UEDBase_H */
