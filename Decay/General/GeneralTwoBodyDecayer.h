// -*- C++ -*-
//
// GeneralTwoBodyDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralTwoBodyDecayer_H
#define HERWIG_GeneralTwoBodyDecayer_H
//
// This is the declaration of the GeneralTwoBodyDecayer class.
//

#include "Herwig/Decay/PerturbativeDecayer.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "GeneralTwoBodyDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

/** \ingroup Decay
 * The GeneralTwoBodyDecayer class is designed to be the base class 
 * for 2 body decays for some general model. It inherits from 
 * PerturbativeDecayer and implements the modeNumber() virtual function
 * that is the  same for all of the decays. A decayer for
 * a specific spin configuration should inherit from this and implement
 * the me2() and partialWidth() member functions. The colourConnections()
 * member should be called from inside me2() in the inheriting decayer
 * to set up the colour lines.
 *
 * @see \ref GeneralTwoBodyDecayerInterfaces "The interfaces"
 * defined for GeneralTwoBodyDecayer.
 * @see PerturbativeDecayer
 */
class GeneralTwoBodyDecayer: public PerturbativeDecayer {

public:
  
  /** A ParticleData ptr and (possible) mass pair.*/
  typedef pair<tcPDPtr, Energy> PMPair;

public:

  /**
   * The default constructor.
   */
  GeneralTwoBodyDecayer() : maxWeight_(1.), colour_(1,DVector(1,1.))
  {}


  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
			       const tPDVector & children) const;

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent,const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int , const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const = 0;
  
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

  /**
   *  Set the information on the decay
   */
  virtual void setDecayInfo(PDPtr incoming, PDPair outgoing,
			    vector<VertexBasePtr>,
			    map<ShowerInteraction,VertexBasePtr> &,
			    const vector<map<ShowerInteraction,VertexBasePtr> > &,
			    map<ShowerInteraction,VertexBasePtr>) =0;

protected:
  
  /** @name Functions used by inheriting decayers. */
  //@{
  /**
   * Set integration weight
   * @param wgt Maximum integration weight 
   */
  void setWeight(double wgt) { maxWeight_ = wgt; }

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

  /**
   *  Calculate matrix element ratio R/B
   */
  double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
			    const ParticleVector & decay3, MEOption meopt,
			    ShowerInteraction inter);

  /**
   *  Set the information on the decay
   */
  void decayInfo(PDPtr incoming, PDPair outgoing);
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
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  Member for the generation of additional hard radiation
   */
  //@{
  /**
   * Return the matrix of colour factors 
   */
  typedef vector<pair<int,double > > CFlowPairVec;
  typedef vector<CFlowPairVec> CFlow;

  const vector<DVector> & getColourFactors(const Particle & inpart, 
					   const ParticleVector & decay, 
					   unsigned int & nflow); 
 
  const CFlow & colourFlows(const Particle & inpart,
			    const ParticleVector & decay);

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay,
			     ShowerInteraction inter, MEOption meopt);

  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralTwoBodyDecayer & operator=(const GeneralTwoBodyDecayer &) = delete;
 
private:

  /**
   *  Store the incoming particle
   */
  PDPtr incoming_;

  /**
   *  Outgoing particles
   */
  vector<PDPtr> outgoing_;

  /**
   * Maximum weight for integration
   */
  double maxWeight_;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> colour_;

};


/**
 * Write a map with ShowerInteraction as the key
 */
template<typename T, typename Cmp, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const map<ShowerInteraction,T,Cmp,A> & m) {
  os << m.size();
  if(m.find(ShowerInteraction::QCD)!=m.end()) {
    os << 0 << m.at(ShowerInteraction::QCD);
  }
  if(m.find(ShowerInteraction::QED)!=m.end()) {
    os << 1 << m.at(ShowerInteraction::QED);
  }
  return os;
}

/**
 * Read a map with ShowerInteraction as the key
 */
template <typename T, typename Cmp, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, map<ShowerInteraction,T,Cmp,A> & m) {
  m.clear();
  long size;
  int k;
  is >> size;
  while ( size-- && is ) {
    is >> k;
    if(k==0)
      is >> m[ShowerInteraction::QCD];
    else if(k==1)
      is >> m[ShowerInteraction::QED];
    else
      assert(false);
  }
  return is;
}
}

#endif /* HERWIG_GeneralTwoBodyDecayer_H */
