// -*- C++ -*-
#ifndef HERWIG_EvtGen_H
#define HERWIG_EvtGen_H
//
// This is the declaration of the EvtGen class.
//

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtRaritaSchwinger.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenRandom.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Helicity/LorentzSpinor.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/RhoDMatrix.h"
#include "ThePEG/Helicity/LorentzRSSpinor.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/EventRecord/Particle.h"
#include "EvtGen.fh"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/Helicity/TensorSpinInfo.h"
#include "ThePEG/Helicity/SpinInfo.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The EvtGen class is the main class for the use of the EvtGen decay
 * package with Herwig++. It is designed to replace the EvtGen class of the 
 * EvtGen package and provide additional functionality needed by Herwig++.
 *
 * @see \ref EvtGenInterfaces "The interfaces"
 * defined for EvtGen.
 */
class EvtGen: public Interfaced {

public:

  /**
   * The default constructor.
   */
  inline EvtGen();

  /**
   *  Use EvtGen to perform a decay
   * @param parent The decaying particle
   * @param recursive Whether or not EvtGen should recursively decay the
   * products of the decay
   * @param dm The decaymode
   * @return The decay products
   */
  ParticleVector decay(const Particle &parent,bool recursive,
		       const DecayMode & dm) const;
  
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
  static ClassDescription<EvtGen> initEvtGen;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EvtGen & operator=(const EvtGen &);

private:

  /** @name Functions to convert between EvtGen and Herwig++ classes */
  //@{
  /**
   * Convert a particle to an EvtGen particle.
   * @param part The particle to be converted.
   */
  EvtParticle *EvtGenParticle(const Particle & part) const;

  /**
   * Convert a particle from an EvtGen one to ThePEG one.
   * @param evtpart The EvtGen particle.
   * @param pd Pointer to the particle data object of ThePEG for the particle.
   * @param spin Convert the spin information as well
   */
  inline PPtr ThePEGParticle(EvtParticle *evtpart, tcPDPtr pd,bool spin) const;

  /**
   * Check the particle has SpinInfo and if not create it
   * @param part The particle
   */
  inline tSpinfoPtr getSpinInfo(const Particle &part) const;

  /**
   * Set the SpinInfo for a ThePEG particle using an EvtGen particle
   * @param pegpart ThePEG particle.
   * @param evtpart The EvtGen particle.
   */
  void ThePEGSpin(PPtr pegpart,EvtParticle *evtpart) const;

  /**
   *  Return the decay products of an EvtGen particle in as ThePEG particles
   * @param evtpart The EvtGen particle
   * @param spin Produce spin information for the particles
   */
  ParticleVector decayProducts(EvtParticle* evtpart,bool spin) const;

  /**
   * Convert a Lorentz5Momentum to a real EvtGen 4-vector
   * @param mom The momentum to be converted
   */
  inline EvtVector4R EvtGenMomentum(const Lorentz5Momentum & mom) const;

  /**
   * Convert from EvtGen momentum to Lorentz5Momentum
   * @param mom The EvtGen 4-momentum
   * @param mass The mass
   */
  inline Lorentz5Momentum ThePEGMomentum(const EvtVector4R & mom,double mass) const;

  /**
   * Convert a spin density matrix to an EvtGen spin density matrix.
   * @param rho The spin density matrix to be converted.
   */
  inline EvtSpinDensity EvtGenSpinDensity(const RhoDMatrix & rho) const;

  /**
   * Convert a spin density to a ThePEG one from an EvtGen one
   * @param rho The spin density matrix to be converted
   * @param id The PDG code of the particle to get special cases right.
   */
  RhoDMatrix ThePEGSpinDensity(const EvtSpinDensity & rho, int id) const;

  /**
   * Convert from our complex to the EvtGen one
   */
  inline EvtComplex EvtGenComplex(Complex) const;

  /**
   * Convert from EvtGen complex to ours
   */
  inline Complex ThePEGComplex(EvtComplex) const;

  /**
   * Convert a LorentzSpinor to an EvtGen one. The spinor is converted to the 
   * EvtGen Dirac representation/
   * @param sp The LorentzSpinor
   */
  inline EvtDiracSpinor EvtGenSpinor(const LorentzSpinor<SqrtEnergy> & sp) const;

  /**
   * Convert an EvtDiracSpinor a LorentzSpinor. This spinor is converted to 
   * the default Dirac matrix representation used by ThePEG.
   * @param sp The EvtDiracSpinor
   */
  inline LorentzSpinor<SqrtEnergy> ThePEGSpinor(const EvtDiracSpinor & sp) const;

  /**
   * Convert a LorentzPolarizationVector to a complex EvtGen 4-vector
   * @param eps The polarization vector to be converted
   */
  inline EvtVector4C EvtGenPolarization(const LorentzPolarizationVector & eps) const;

  /**
   * Convert an EvtGen complex 4-vector to a LorentzPolarizationVector
   * @param eps The complex 4-vector to be converted.
   */
  inline LorentzPolarizationVector ThePEGPolarization(const EvtVector4C & eps) const;

  /**
   * Convert our Rarita-Schwinger spinor to the EvtGen one
   * @param sp Our  RS Spinor
   */
  inline EvtRaritaSchwinger EvtGenRSSpinor(const LorentzRSSpinor<SqrtEnergy> & sp) const;

  /**
   * Convert an EvtGen Rarita-Schwinger spinor to ours
   * @param sp The EvtGen RS spinor.
   */
  inline LorentzRSSpinor<SqrtEnergy> ThePEGRSSpinor(const EvtRaritaSchwinger & sp) const;

  /**
   * Convert our tensor to the EvtGen one.
   * @param ten Our tensor
   */
  inline EvtTensor4C EvtGenTensor(const LorentzTensor<double> & ten) const;

  /**
   * Convert an EvtGen tensor to ThePEG
   * @param ten The EvtGen tensor
   */
  inline LorentzTensor<double> ThePEGTensor(const EvtTensor4C & ten) const;

  /**
   * Convert a PDG code from ThePEG into an EvtGen particle id
   * @param id The PDG code
   * @param exception Whether or not to throw an Exception if fails
   */
  EvtId EvtGenID(int id,bool exception=true) const;

  /**
   * Convert an EvtGen EvtId to a PDG code in our conventions
   * @param id The EvtGen ID.
   * @param exception Whether or not to throw an Exception if fails
   */
  int   ThePEGID(EvtId id,bool exception=true) const;

  /**
   *  Construct the DecayVertex for Herwig using the information from
   *  EvtGen
   * @param parent The decaying particle
   * @param products The outgoing particles
   * @param damp Pointer to the EvtGen decayer
   */
  void constructVertex(const Particle & parent,ParticleVector products,
		       EvtDecayAmp* damp) const;
  //@}

  /**
   *  Find the location in the EvtGen list of decay channels for
   *  a given decay mode.
   */
  int EvtGenChannel(const DecayMode &dm) const;

  /**
   *  Check the conversion of particles between Herwig++ and EvtGen
   */
  void checkConversion() const;

  /**
   * Output the EvtGen decay modes for a given particle
   * @param id The PDG code of the particle to output
   */
  void outputEvtGenDecays(long id) const;

private:

  /**
   *  The name of the file containing the decays
   */
  string _decayname;

  /**
   *  The name of the file containing the particle data
   */
  string _pdtname;

  /**
   *  The EvtGen particle data tables
   */
  EvtPDL _pdl;

  /**
   * Pointer to the random number generator for EvtGen
   */
  EvtRandomEngine *_evtrnd;

  /**
   *  Maximum number of attempts for EvtGen to generate a decay
   */
  unsigned int _maxtry;

  /**
   *  Maximum number of tries for EvtGen to unweight
   */
  unsigned int _maxunwgt;

  /**
   *  Check the conversion of the particles
   */
  bool _checkconv;

  /**
   *  Particles for which to output the EvtGen decays
   */
  vector<long> _convid;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EvtGen. */
template <>
struct BaseClassTrait<Herwig::EvtGen,1> {
  /** Typedef of the first base class of EvtGen. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EvtGen class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EvtGen>
  : public ClassTraitsBase<Herwig::EvtGen> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EvtGen"; }
  /** Return the name of the shared library be loaded to get
   *  access to the EvtGen class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwEvtGen.so"; }
};

/** @endcond */

}

#include "EvtGen.icc"

#endif /* HERWIG_EvtGen_H */
