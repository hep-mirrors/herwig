// -*- C++ -*-
#ifndef HERWIG_GeneralThreeBodyDecayer_H
#define HERWIG_GeneralThreeBodyDecayer_H
//
// This is the declaration of the GeneralThreeBodyDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Models/General/TBDiagram.h"
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
  inline GeneralThreeBodyDecayer();

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
  void setDecayInfo(PDPtr incoming,vector<PDPtr> outgoing,
		    const vector<TBDiagram> & process,
		    const vector<DVector> & factors,
		    const vector<DVector> & Ncfactors,
		    const unsigned int ncf);

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
  //@}

protected:

  /**
   * Access the TBDiagrams that store the required information
   * to create the diagrams
   */
  inline const vector<TBDiagram> & getProcessInfo() const;

  /**
   *  Incoming particle
   */
  inline PDPtr incoming() const;

  /**
   *  Outgoing particles
   */
  inline const vector<PDPtr> & outgoing() const;
 
  /**
   *  Number of colour flows
   */
  inline unsigned int numberOfFlows() const;

  /**
   * Return the matrix of colour factors 
   */
  inline const vector<DVector> & getColourFactors() const;

  /**
   * Return the matrix of colour factors 
   */
  inline const vector<DVector> & getLargeNcColourFactors() const;

  /**
   *  Get the mapping between the phase-space channel and the diagram
   */
  inline const vector<unsigned int> & diagramMap() const;

  /**
   *  Option for the handling of the widths of the intermediate particles
   */
  inline unsigned int widthOption() const;

  /**
   * Set colour connections
   * @param parent Parent particle
   * @param out Particle vector containing particles to 
   * connect colour lines
   */
  void colourConnections(const Particle & parent, 
			 const ParticleVector & out) const;

  /**
   *
   */
  void constructIntegratorChannels(vector<int> & intype, vector<Energy> & inmass,
				   vector<Energy> & inwidth, vector<double> & inpow,
				   vector<double> & inweights) const;

  /**
   *  Set the colour flow
   */
  inline void colourFlow(unsigned int) const;

  /**
   *  Set the colour flow
   */
  inline unsigned int const & colourFlow() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<GeneralThreeBodyDecayer> initGeneralThreeBodyDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralThreeBodyDecayer & operator=(const GeneralThreeBodyDecayer &);

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
  vector<TBDiagram> _diagrams;

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
   *  Reference to object to calculate the partial width
   */
  mutable WidthCalculatorBasePtr _widthcalc;

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

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralThreeBodyDecayer. */
template <>
struct BaseClassTrait<Herwig::GeneralThreeBodyDecayer,1> {
  /** Typedef of the first base class of GeneralThreeBodyDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralThreeBodyDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralThreeBodyDecayer>
  : public ClassTraitsBase<Herwig::GeneralThreeBodyDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GeneralThreeBodyDecayer"; }
};

/** @endcond */

}

#include "GeneralThreeBodyDecayer.icc"

#endif /* HERWIG_GeneralThreeBodyDecayer_H */
