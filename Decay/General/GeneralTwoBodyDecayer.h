// -*- C++ -*-
//
// GeneralTwoBodyDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralTwoBodyDecayer_H
#define HERWIG_GeneralTwoBodyDecayer_H
//
// This is the declaration of the GeneralTwoBodyDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "GeneralTwoBodyDecayer.fh"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

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
  GeneralTwoBodyDecayer() : _maxweight(1.), mb_(ZERO), e_(0.), s_(0.), e2_(0.), s2_(0.), 
			    pTmin_(GeV), pT_(ZERO), colour_(1,DVector(1,1.))
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

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}

  /**
   *  Member to generate the hardest emission in the POWHEG scheme
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);


  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay, MEOption meopt);
  //@}

  /**
   *  Set the information on the decay
   */
  void setDecayInfo(PDPtr incoming,PDPair outgoing,
		    VertexBasePtr,VertexBasePtr,
		    const vector<VertexBasePtr> &,
		    VertexBasePtr);

protected:
  
  /** @name Functions used by inheriting decayers. */
  //@{
  
  /**
   * Get vertex pointer
   * @return a pointer to the vertex
   */
  VertexBasePtr getVertex() const { return vertex_; }

  /**
   * Get vertex pointer
   * @return a pointer to the vertex for QCD radiation off the decaying particle
   */
  VertexBasePtr getIncomingVertex() const { return  incomingVertex_; }

  /**
   * Get vertex pointer
   * @return a pointer to the vertex for QCD radiation off the decay products
   */
  vector<VertexBasePtr> getOutgoingVertices() const { return outgoingVertices_; }

  /**
   * Get vertex pointer
   * @return a pointer to the vertex for QCD radiation from 4 point vertex
   */
  VertexBasePtr getFourPointVertex() const { return fourPointVertex_; }


  /**
   * Set integration weight
   * @param wgt Maximum integration weight 
   */
  void setWeight(double wgt) { _maxweight = wgt; }

  /**
   * Set colour connections
   * @param parent Parent particle
   * @param out Particle vector containing particles to 
   * connect colour lines
   */
  void colourConnections(const Particle & parent, 
			 const ParticleVector & out) const;

  /**
   * Type of dipole
   */
  enum dipoleType {FFa, FFc, IFa, IFc, IFba, IFbc};

  /**
   *  Compute the spin and colour factor
   */
  double colourFactor(tcPDPtr in, tcPDPtr out1, tcPDPtr out2) const;

  /**
   *  Calculate matrix element ratio R/B
   */
  double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
			    const ParticleVector & decay3, MEOption meopt);

  /**
   *  Calculate momenta of all the particles
   */
  bool calcMomenta(int j, Energy pT, double y, double phi, double& xg, 
		   double& xs, double& xe, double& xe_z, 
		   vector<Lorentz5Momentum>& particleMomenta);

  /**
   *  Check the calculated momenta are physical
   */
  bool psCheck(const double xg, const double xs);

  /**
   *  Return the momenta including the hard emission
   */
  vector<Lorentz5Momentum> hardMomenta(const PPtr &in, 
				       const PPtr &emitter, 
				       const PPtr &spectator, 
				       const vector<dipoleType>  &dipoles, int i);

  /**
   * Return dipole corresponding to the dipoleType dipoleId
   */
  InvEnergy2 calculateDipole(const dipoleType & dipoleId,   const Particle & inpart,
			     const ParticleVector & decay3, const dipoleType & emittingDipole);

  /**
   * Return contribution to dipole that depends on the spin of the emitter
   */
  double dipoleSpinFactor(const PPtr & emitter, double z);

  /**
   *  Work out the type of process
   */
  bool identifyDipoles(vector<dipoleType> & dipoles,
		       PPtr & aProgenitor,
		       PPtr & bProgenitor,
		       PPtr & cProgenitor) const;
  
  /**
   * Set up the colour lines
   */
  void getColourLines(RealEmissionProcessPtr real);


  /**
   *  Return the colour coefficient of the dipole
   */
  double colourCoeff(const PDT::Colour emitter, const PDT::Colour spectator,
		     const PDT::Colour other);

  /**
   *  Coupling for the generation of hard radiation
   */
  ShowerAlphaPtr coupling() {return coupling_;}
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
   *  Store the incoming particle
   */
  PDPtr _incoming;

  /**
   *  Outgoing particles
   */
  vector<PDPtr> _outgoing;
  
  /**
   * Pointer to vertex
   */
  VertexBasePtr vertex_;

  /**
   *  Pointer to vertex for radiation from the incoming particle
   */
  VertexBasePtr incomingVertex_;

  /**
   *  Pointer to the vertices for radiation from the outgoing particles
   */
  vector<VertexBasePtr> outgoingVertices_; 

  /**
   *  Pointer to vertex for radiation coming from 4 point vertex
   */
  VertexBasePtr fourPointVertex_;


  /**
   * Maximum weight for integration
   */
  double _maxweight;

  /**
   *  Mass of decaying particle
   */
  Energy mb_;

  /**
   *  Reduced mass of emitter child particle
   */
  double e_;

  /**
   * Reduced mass of spectator child particle
   */
  double s_;

  /**
   *  Reduced mass of emitter child particle squared
   */
  double e2_;

  /**
   * Reduced mass of spectator child particle squared
   */
  double s2_;

  /**
   *  Minimum \f$p_T\f$
   */
  Energy pTmin_;

  /**
   *  Transverse momentum of the emission
   */
  Energy pT_;

  /**
   *  Coupling for the generation of hard radiation
   */
  ShowerAlphaPtr coupling_;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> colour_;

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
