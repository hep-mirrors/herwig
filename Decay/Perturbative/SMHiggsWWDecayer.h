// -*- C++ -*-
#ifndef HERWIG_SMHiggsWWDecayer_H
#define HERWIG_SMHiggsWWDecayer_H
//
// This is the declaration of the SMHiggsWWDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.fh"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.fh"
#include "SMHiggsWWDecayer.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The SMHiggsWWDecayer class performs the decay of the Standard Model
 * Higgs boson to \f$W^+W^-\f$ and \f$Z^0Z^0\f$ including the decays
 * of the gauge bosons.
 *
 * @see \ref SMHiggsWWDecayerInterfaces "The interfaces"
 * defined for SMHiggsWWDecayer.
 */
class SMHiggsWWDecayer: public DecayIntegrator {

public:

  /**
   *  A typedef to select the  boson decay modes
   */
  typedef Selector<unsigned int> ModeSelector;

public:

  /**
   * The default constructor.
   */
  SMHiggsWWDecayer();

  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const PDVector & children) const;

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
   */
  virtual int modeNumber(bool &, tcPDPtr, const PDVector & ) const {return -1;}
  
  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMHiggsWWDecayer> initSMHiggsWWDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsWWDecayer & operator=(const SMHiggsWWDecayer &);

private:

  /**
   *  Pointers to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the fermion-femion-W vertex
   */
  FFVVertexPtr _theFFWVertex;

  /**
   *  Pointer to the fermion-femion-Z vertex
   */
  FFVVertexPtr _theFFZVertex;

  /**
   *  Pointer to the higgs-WW/ZZ vertex
   */
  VVSVertexPtr _theHVVVertex;
  //@}

  /**
   *  Selectors for the gauge boson decay modes
   */
  //@{
  /**
   *  Selector for the W decays
   */
  ModeSelector _wdecays;

  /**
   *  Selector for the Z decays
   */
  ModeSelector _zdecays;
  //@}

  /**
   *  Product of gauge boson branching ratios for normalisation
   */
  vector<double> _ratio;

  /**
   *  Maximum weights for the decays
   */
  //@{
  /**
   *  Maximum weight for \f$H\to W^+W^-\f$ decays
   */
  double _wmax;

  /**
   *  Maximum weight for \f$H\to Z^0Z^0\f$ decays
   */
  double _zmax;
  //@}

  /**
   *  Parameters for the generation of the boson masses
   */
  //@{
  /**
   *  Generate using a Breit-Wigner or a power law
   */
  bool _breit;

  /**
   *  Power for the power law
   */
  double _power;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsWWDecayer. */
template <>
struct BaseClassTrait<Herwig::SMHiggsWWDecayer,1> {
  /** Typedef of the first base class of SMHiggsWWDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsWWDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsWWDecayer>
  : public ClassTraitsBase<Herwig::SMHiggsWWDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsWWDecayer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SMHiggsWWDecayer is implemented. It may also include several, space-separated,
   * libraries if the class SMHiggsWWDecayer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPerturbativeHiggsDecay.so"; }
};

/** @endcond */

}

#include "SMHiggsWWDecayer.icc"

#endif /* HERWIG_SMHiggsWWDecayer_H */
