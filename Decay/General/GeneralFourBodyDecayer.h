// -*- C++ -*-
#ifndef HERWIG_GeneralFourBodyDecayer_H
#define HERWIG_GeneralFourBodyDecayer_H
//
// This is the declaration of the GeneralFourBodyDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Models/General/PrototypeVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GeneralFourBodyDecayer class.
 *
 * @see \ref GeneralFourBodyDecayerInterfaces "The interfaces"
 * defined for GeneralFourBodyDecayer.
 */
class GeneralFourBodyDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  GeneralFourBodyDecayer(): _nflow(999), _widthopt(1), 
			    _reftag(), _reftagcc(), _iflow(999)
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
  virtual int modeNumber(bool & cc, tcPDPtr parent,
			 const tPDVector & children) const;

  /**
   *  Set the diagrams
   */
  bool setDecayInfo(PDPtr incoming,vector<PDPtr> outgoing,
		    const vector<NBDiagram> & process,
		    double symfac);
  //@}
  
  /**
   * Function to return partial Width
   * @param inpart Pointer to incoming particle data object
   * @param outgoing the decay products
   */
  virtual Energy partialWidth(tPDPtr inpart,
			      OrderedParticles outgoing) const;

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
   *  Incoming particle
   */
  PDPtr incoming() const { return _incoming; }
  
  /**
   *  Outgoing particles
   */
  const vector<PDPtr> & outgoing() const {  return _outgoing; }
  
  /**
   *  Number of colour flows
   */
  unsigned int numberOfFlows() const { return _nflow; }
  
  /**
   * Set up the colour factors
   */
  bool setColourFactors(double symfac);
  
  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getColourFactors() const {  return _colour; }
  
  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getLargeNcColourFactors() const {
    return _colourLargeNC;
  }
  
  /**
   *  Option for the handling of the widths of the intermediate particles
   */
  unsigned int widthOption() const { return _widthopt; }
  
  /**
   * Set colour connections
   * @param parent Parent particle
   * @param out Particle vector containing particles to 
   * connect colour lines
   */
  void colourConnections(const Particle & parent, 
			 const ParticleVector & out) const;

  /**
   *  Set the colour flow
   * @param flow The value for the colour flow
   */
  void colourFlow(unsigned int flow) const { _iflow = flow; }
  
  /**
   *  Set the colour flow
   */
  unsigned int const & colourFlow() const { return _iflow; }

  /**
   * Access the TBDiagrams that store the required information
   * to create the diagrams
   */
  const vector<NBDiagram> & getProcessInfo() const {
    return _diagrams;
  }

  /**
   *  Get the mapping between the phase-space channel and the diagram
   */
  const vector<unsigned int> & diagramMap() const { 
    return _diagmap; 
  }
  
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
  GeneralFourBodyDecayer & operator=(const GeneralFourBodyDecayer &);

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
   *  Store the diagrams for the decay
   */
  vector<NBDiagram> _diagrams;

  /**
   *  Map between the diagrams and the phase-space channels
   */
  vector<unsigned int> _diagmap;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> _colour;

  /**
   *  Store cololur factors for ME calc at large N_c
   */
  vector<DVector> _colourLargeNC;

  /**
   * The number of colourflows.
   */
  unsigned int _nflow;

  /**
   *  Option for the treatment of the widths 
   */
  unsigned int _widthopt;

  /**
   * Store a decay tag for this mode that can be tested when
   * trying to determine whether it can be generated by
   * this Decayer
   */
  string _reftag;

  /**
   * Store a decay tag for the cc-mode that can be tested when
   * trying to determine whether it can be generated by
   * this Decayer
   */
  string _reftagcc;

  /**
   *  The colour flow
   */
  mutable unsigned int _iflow;
};

}

#endif /* HERWIG_GeneralFourBodyDecayer_H */
