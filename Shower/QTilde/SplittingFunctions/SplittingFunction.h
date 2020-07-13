// -*- C++ -*-
//
// SplittingFunction.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SplittingFunction_H
#define HERWIG_SplittingFunction_H
//
// This is the declaration of the SplittingFunction class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Shower/QTilde/ShowerConfig.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "Herwig/Decay/DecayMatrixElement.h"
#include "Herwig/Shower/QTilde/Kinematics/ShowerKinematics.fh"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/PDT/ParticleData.h"
#include "SplittingFunction.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This is an abstract class which defines the common interface
 *  for all \f$1\to2\f$ splitting functions, for both initial-state
 *  and final-state radiation. 
 *
 *  The SplittingFunction class contains a number of purely virtual members
 *  which must be implemented in the inheriting classes. The class also stores
 *  the interaction type of the spltting function.
 *
 *  The inheriting classes need to specific the splitting function 
 *  \f$P(z,2p_j\cdot p_k)\f$, in terms of the energy fraction \f$z\f$ and
 *  the evolution scale. In order to allow the splitting functions to be used
 *  with different choices of evolution functions the scale is given by
 * \f[2p_j\cdot p_k=(p_j+p_k)^2-m_{jk}^2=Q^2-(p_j+p_k)^2=z(1-z)\tilde{q}^2=
 *   \frac{p_T^2}{z(1-z)}-m_{jk}^2+\frac{m_j^2}{z}+\frac{m_k^2}{1-z},\f]
 *  where \f$Q^2\f$ is the virtuality of the branching particle,
 *  $p_T$ is the relative transverse momentum of the branching products and
 *  \f$\tilde{q}^2\f$ is the angular variable described in hep-ph/0310083.
 *
 *  In addition an overestimate of the 
 *  splitting function, \f$P_{\rm over}(z)\f$ which only depends upon \f$z\f$, 
 *  the integral, inverse of the integral for this overestimate and
 *  ratio of the true splitting function to the overestimate must be provided
 *  as they are necessary for the veto alogrithm used to implement the evolution.
 *
 * @see \ref SplittingFunctionInterfaces "The interfaces"
 * defined for SplittingFunction.
 */
class SplittingFunction: public Interfaced {
  
  friend class SudakovFormFactor;
  
public:

  /**
   * The default constructor.   
   * @param b All splitting functions must have an interaction order
   */
  SplittingFunction()
    : Interfaced() {}

public:

  /**
   *  Methods to return the interaction type and order for the splitting function
   */
  //@{
  /**
   *  Purely virtual method which should determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const = 0;
  //@}

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * Purely virtual method which should return the exact value of the splitting function,
   * \f$P\f$ evaluated in terms of the energy fraction, \f$z\f$, and the evolution scale 
   \f$\tilde{q}^2\f$.
   * @param z   The energy fraction.
   * @param t   The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids,
		   const bool mass, const RhoDMatrix & rho) const = 0;

  /**
   * Purely virtual method which should return
   * an overestimate of the splitting function,
   * \f$P_{\rm over}\f$ such that the result \f$P_{\rm over}\geq P\f$. This function
   * should be simple enough that it does not depend on the evolution scale.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList & ids) const = 0; 

  /**
   * Purely virtual method which should return
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			const bool mass, const RhoDMatrix & rho) const = 0;

  /**
   * Purely virtual method which should return the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z         The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   *                  
   */
  virtual double integOverP(const double z, const IdList & ids, 
			    unsigned int PDFfactor=0) const = 0; 

  /**
   * Purely virtual method which should return the inverse of the 
   * indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$ which is used to
   * generate the value of \f$z\f$.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double invIntegOverP(const double r, const IdList & ids, 
			       unsigned int PDFfactor=0) const = 0;
  //@}

  /**
   * Method to calculate the azimuthal angle for forward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiForward(const double z, const Energy2 t, const IdList & ids,
		     const RhoDMatrix &) = 0;

  /**
   * Method to calculate the azimuthal angle for backward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiBackward(const double z, const Energy2 t, const IdList & ids,
		      const RhoDMatrix &) = 0;

  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param phi The azimuthal angle, \f$\phi\f$.
   * @param timeLike Whether timelike or spacelike, affects inclusive of mass terms
   */
  virtual DecayMEPtr matrixElement(const double z, const Energy2 t, 
				   const IdList & ids, const double phi,
                                   bool timeLike) = 0;

public:

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
  SplittingFunction & operator=(const SplittingFunction &) = delete;

};

}

#endif /* HERWIG_SplittingFunction_H */
