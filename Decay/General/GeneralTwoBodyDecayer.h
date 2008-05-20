// -*- C++ -*-
//
// GeneralTwoBodyDecayer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralTwoBodyDecayer_H
#define HERWIG_GeneralTwoBodyDecayer_H
//
// This is the declaration of the GeneralTwoBodyDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "GeneralTwoBodyDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

/** \ingroup Decay
 * The GeneralTwoBodyDecayer class is designed to be the base class 
 * for 2 body decays for some general model. It inherits from 
 * DecayIntegrator and implements the modeNumber() virtual function
 * that is the  same for all of the decays. A decayer for
 * a specific spin configuration should inherit from this and implement
 * the me2() and partialWidth() member functions. The colourConnections()
 * member should be called from inside me2() in the inheriting decayer
 * to set up the colour lines.
 *
 * @see \ref GeneralTwoBodyDecayerInterfaces "The interfaces"
 * defined for GeneralTwoBodyDecayer.
 * @see DecayIntegrator
 */
class GeneralTwoBodyDecayer: public DecayIntegrator {

public:
  
  /** A ParticleData ptr and (possible) mass pair.*/
  typedef pair<tcPDPtr, Energy> PMPair;

public:

  /**
   * The default constructor.
   */
  inline GeneralTwoBodyDecayer() : _thelist(0,0), _maxweight(1,1.) {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const PDVector & children) const;

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent,const PDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const = 0;
  
  /**
   * Function to return partial Width
   * @param inpart The decaying particle.
   * @param outa One of the decay products.
   * @param outb The other decay product.
   */
  virtual Energy partialWidth(PMPair inpart, PMPair outa, 
			      PMPair outb) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width 
   * calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class.
   * @param coupling The coupling for the matrix element.
   * @return True if the the order of the particles in the 
   * decayer is the same as the DecayMode tag.
   */
  virtual bool twoBodyMEcode(const DecayMode & dm, int & mecode,
			     double & coupling) const;

  /**
   * An overidden member to calculate a branching ratio for a certain
   * particle instance.
   * @param dm The DecayMode of the particle
   * @param p The particle object
   * @param oldbrat The branching fraction given in the DecayMode object
   */
  virtual double brat(const DecayMode & dm, const Particle & p,
		      double oldbrat) const;
  //@}

protected:
  
  /** @name Functions used by inheriting decayers. */
  //@{
  /** Set list to search
   * @param ilist 
   */
  inline void addToSearchList(unsigned int ilist) { 
    _thelist.push_back(ilist); 
  }
  
  /**
   * Get vertex pointer
   * @return a pointer to the vertex
   */
  inline VertexBasePtr getVertex() const { return _theVertex; }

  /**
   * Set integration weight
   * @param wgt Maximum integration weight 
   */
  inline void setWeight(const vector<double> & wgt) { _maxweight = wgt; }

  /**
   * Set colour connections
   * @param parent Parent particle
   * @param out Particle vector containing particles to 
   * connect colour lines
   */
  void colourConnections(const Particle & parent, 
			 const ParticleVector & out) const;

  /**
   *  Compute the spin and colour factor
   */
  double colourFactor(tcPDPtr in, tcPDPtr out1, tcPDPtr out2) const;
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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<GeneralTwoBodyDecayer> initGeneralTwoBodyDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralTwoBodyDecayer & operator=(const GeneralTwoBodyDecayer &);
 
private:

  /**
   * vector of ints as to which list(s) to search
   */
  vector<unsigned int> _thelist;
  
  /**
   * Pointer to vertex set in inheriting class
   */
  VertexBasePtr _theVertex;
  
  /**
   * PDG codes for all incoming particles
   **/
  vector<int> _inpart;

  /**
   * PDG codes for 1st set of outgoing particles
   **/
  vector<int> _outparta;
   
  /**
   * PDG codes for 2nd set of outgoing particles
   **/
  vector<int> _outpartb;
 
  /**
   * Vector of maximum weights for integration
   */
  vector<double> _maxweight;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralTwoBodyDecayer. */
template <>
struct BaseClassTrait<Herwig::GeneralTwoBodyDecayer,1> {
  /** Typedef of the first base class of GeneralTwoBodyDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralTwoBodyDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralTwoBodyDecayer>
  : public ClassTraitsBase<Herwig::GeneralTwoBodyDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GeneralTwoBodyDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_GeneralTwoBodyDecayer_H */
