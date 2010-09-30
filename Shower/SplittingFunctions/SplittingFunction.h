// -*- C++ -*-
//
// SplittingFunction.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SplittingFunction_H
#define HERWIG_SplittingFunction_H
//
// This is the declaration of the SplittingFunction class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "Herwig++/Decay/DecayMatrixElement.h"
#include "SplittingFunction.fh"

namespace Herwig {

using namespace ThePEG;

  /** \ingroup Shower
   * Enum to define the possible types of colour structure which can occur in
   * the branching.
   */
  enum ColourStructure {Undefined=-1,TripletTripletOctet=0,OctetOctetOctet=1,
			OctetTripletTriplet=2,TripletOctetTriplet=3};


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

public:

  /**
   * The default constructor.   
   * @param b All splitting functions must have an interaction order
   */
  SplittingFunction(unsigned int b)
    : Interfaced(), _interactionType(ShowerInteraction::UNDEFINED),
      _interactionorder(b), 
      _colourStructure(Undefined), _colourFactor(-1.) {}

public:

  /**
   *  Methods to return the interaction type and order for the splitting function
   */
  //@{
  /**
   *  Return the type of the interaction
   */
  ShowerInteraction::Type interactionType() const {return _interactionType;}

  /**
   *  Return the order of the splitting function in the interaction
   */
  unsigned int interactionOrder() const {return _interactionorder;}

  /**
   *  Return the colour structure
   */
  ColourStructure colourStructure() const {return _colourStructure;}

  /**
   *  Return the colour factor
   */
  double colourFactor() const {return _colourFactor;}
  //@}

  /**
   *  Purely virtual method which should determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const = 0;

  /**
   *  Method to check the colours are correct
   */
  virtual bool checkColours(const IdList & ids) const;

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
		   const bool mass) const = 0;

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
			const bool mass) const = 0;

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
   * Purely virtual method which should make the proper colour connection 
   * between the emitting parent and the branching products.
   * @param parent The parent for the branching
   * @param first  The first  branching product
   * @param second The second branching product
   * @param back Whether this is foward or backward evolution.
   */
  virtual void colourConnection(tShowerParticlePtr parent,
				tShowerParticlePtr first,
				tShowerParticlePtr second,
				const bool back) const;

  /**
   * Method to calculate the azimuthal angle
   * @param particle The particle which is branching
   * @param showerkin The ShowerKinematics object
   * @param z The energy fraction
   * @param qtilde   The scale.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual double generatePhi(ShowerParticle & particle,ShoKinPtr showerkin,
			     const double z, const Energy qtilde, const IdList & ids,
			     const RhoDMatrix &, const double phi);

  /**
   * Calculate the matrix element for the splitting
   */
  virtual DecayMatrixElement matrixElement(ShowerParticle & particle,ShoKinPtr showerkin,
					   const double z, const Energy qtilde, 
					   const IdList & ids,
					   const RhoDMatrix &, const double phi);

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

protected:

  /**
   *  Set the colour factor
   */
  void colourFactor(double in) {_colourFactor=in;}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractClassDescription<SplittingFunction> initSplittingFunction;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SplittingFunction & operator=(const SplittingFunction &);

private:

  /**
   *  The interaction type for the splitting function.
   */
  ShowerInteraction::Type _interactionType;

  /**
   *  The order of the splitting function in the coupling
   */
  unsigned int _interactionorder;

  /**
   *  The colour structure
   */
  ColourStructure _colourStructure;

  /**
   *  The colour factor
   */
  double _colourFactor;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SplittingFunction. */
template <>
struct BaseClassTrait<Herwig::SplittingFunction,1> {
  /** Typedef of the first base class of SplittingFunction. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SplittingFunction class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SplittingFunction>
  : public ClassTraitsBase<Herwig::SplittingFunction> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SplittingFunction"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SplittingFunction is implemented. It may also include several, space-separated,
   * libraries if the class SplittingFunction depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SplittingFunction_H */
