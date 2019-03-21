// -*- C++ -*-
#ifndef HERWIG_GeneralThreeBodyDecayer_H
#define HERWIG_GeneralThreeBodyDecayer_H
//
// This is the declaration of the GeneralThreeBodyDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Models/General/TBDiagram.h"
#include "GeneralThreeBodyDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * Here is the documentation of the GeneralThreeBodyDecayer class.
 *
 * @see \ref GeneralThreeBodyDecayerInterfaces "The interfaces"
 * defined for GeneralThreeBodyDecayer.
 */
class GeneralThreeBodyDecayer: public DecayIntegrator {

public:
  
  /** A ParticleData ptr and (possible) mass pair.*/
  typedef pair<tcPDPtr, Energy> PMPair;


public:

  /**
   * The default constructor.
   */
  GeneralThreeBodyDecayer() : nflow_(999), widthOpt_(1), 
			      refTag_(), refTagCC_(), iflow_(999),
			      intOpt_(0), relerr_(1e-2)
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
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(const int imode,  const Energy2 q2,
					const Energy2 s3, const Energy2 s2, 
					const Energy2 s1, const Energy  m1, 
					const Energy  m2, const Energy  m3) const;
  
  /**
   * Function to return partial Width
   * @param inpart The decaying particle.
   * @param outa First  decay product.
   * @param outb Second decay product.
   * @param outc Third  decay product.
   */
  virtual Energy partialWidth(PMPair inpart, PMPair outa, 
			      PMPair outb, PMPair outc) const;

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
   *  Set the diagrams
   */
  bool setDecayInfo(PDPtr incoming,vector<PDPtr> outgoing,
		    const vector<TBDiagram> & process,
		    double symfac);

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

  /**
   *   Set up the diagrams etc
   */
  virtual void setupDiagrams(bool checkKinematics);

protected:

  /**
   * Access the TBDiagrams that store the required information
   * to create the diagrams
   */
  const vector<TBDiagram> & getProcessInfo() const {
    return diagrams_;
  }

  /**
   *  Incoming particle
   */
  PDPtr incoming() const { return incoming_; }

  /**
   *  Outgoing particles
   */
  const vector<PDPtr> & outgoing() const {  return outgoing_; }
 
  /**
   *  Number of colour flows
   */
  unsigned int numberOfFlows() const { return nflow_; }

  /**
   * Set up the colour factors
   */
  bool setColourFactors(double symfac);

  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getColourFactors() const {  return colour_; }

  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getLargeNcColourFactors() const {
    return colourLargeNC_;
  }

  /**
   *  Get the mapping between the phase-space channel and the diagram
   */
  const vector<unsigned int> & diagramMap() const { 
    return diagmap_; 
  }

  /**
   *  Option for the handling of the widths of the intermediate particles
   */
  unsigned int widthOption() const { return widthOpt_; }

  /**
   * Set colour connections
   * @param parent Parent particle
   * @param out Particle vector containing particles to 
   * connect colour lines
   */
  void colourConnections(const Particle & parent, 
			 const ParticleVector & out) const;

  /**
   *  Method to construct the channels for the integrator to give the partial width
   * @param intype  Types of the channels
   * @param inmass Mass for the channels
   * @param inwidth Width for the channels
   * @param inpow Power for the channels
   * @param inweights Weights for the channels
   */
  void constructIntegratorChannels(vector<int> & intype, vector<Energy> & inmass,
				   vector<Energy> & inwidth, vector<double> & inpow,
				   vector<double> & inweights) const;

  /**
   *  Set the colour flow
   * @param flow The value for the colour flow
   */
  void colourFlow(unsigned int flow) const { iflow_ = flow; }

  /**
   *  Set the colour flow
   */
  unsigned int const & colourFlow() const { return iflow_; }

  /**
   *  Relative error for GQ integration
   */
  double relativeError() const {return relerr_;}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralThreeBodyDecayer & operator=(const GeneralThreeBodyDecayer &) = delete;

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
   *  Store the diagrams for the decay
   */
  vector<TBDiagram> diagrams_;

  /**
   *  Map between the diagrams and the phase-space channels
   */
  vector<unsigned int> diagmap_;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> colour_;

  /**
   *  Store cololur factors for ME calc at large N_c
   */
  vector<DVector> colourLargeNC_;

  /**
   * The number of colourflows.
   */
  unsigned int nflow_;

  /**
   *  Reference to object to calculate the partial width
   */
  mutable WidthCalculatorBasePtr widthCalc_;

  /**
   *  Option for the treatment of the widths 
   */
  unsigned int widthOpt_;

  /**
   * Store a decay tag for this mode that can be tested when
   * trying to determine whether it can be generated by
   * this Decayer
   */
  string refTag_;

  /**
   * Store a decay tag for the cc-mode that can be tested when
   * trying to determine whether it can be generated by
   * this Decayer
   */
  string refTagCC_;

  /**
   *  The colour flow
   */
  mutable unsigned int iflow_;

  /**
   *  Option for the construction of the gaussian integrator
   */
  unsigned int intOpt_;

  /**
   *  Relative error for GQ integration of partial width
   */
  double relerr_;
};

}

#endif /* HERWIG_GeneralThreeBodyDecayer_H */
