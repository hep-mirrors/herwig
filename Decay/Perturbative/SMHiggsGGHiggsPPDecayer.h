// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsGGHiggsPPDecayer_H
#define HERWIG_SMHiggsGGHiggsPPDecayer_H
//
// This is the declaration of the SMHiggsGGHiggsPPDecayer class.
//

#include "Herwig/Decay/PerturbativeDecayer.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
  
/**
 * The <code>SMHiggsGGHiggsPPDecayer</code> class performs the
 * of a Standard Model Higgs boson to:  a pair
 * of photons or a pair of gluons, or a \f$Z^0\f$ boson and a photon.
 *
 * @see PerturbativeDecayer
 */ 
class SMHiggsGGHiggsPPDecayer: public PerturbativeDecayer {
  
public:

  /**
   * The default constructor.
   */
  SMHiggsGGHiggsPPDecayer() : _h0wgt(3,1.), _minloop(6), _maxloop(6), _massopt(0)
  {}
  
  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const;
  
  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   * Which of the possible decays is required
   */
  virtual int modeNumber(bool &, tcPDPtr, const tPDVector & ) const {return -1;}
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,const tPDVector & children) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

  /**
   *  Calculate matrix element ratio R/B
   */
  virtual double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				    const ParticleVector & decay3, MEOption meopt,
				    ShowerInteraction inter);

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}
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
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}
  
  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}
  
protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
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
   *  Calculate the NLO real emission piece of ME
   */
  double realME(const vector<cPDPtr> & partons, 
		const vector<Lorentz5Momentum> & momenta) const;
  
  /**
   *  Calculate the LO ME
   */
  Energy2 loME(Energy mh) const;
  
private:
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsGGHiggsPPDecayer & operator=(const SMHiggsGGHiggsPPDecayer &) = delete;
  
  /**
   * Pointer to h->gluon,gluon vertex
   */
  AbstractVVSVertexPtr _hggvertex;
  
  /**
   * Pointer to h->gamma,gamma vertex
   */
  AbstractVVSVertexPtr _hppvertex;
  
  /**
   * Pointer to h->gamma,gamma vertex
   */
  AbstractVVSVertexPtr _hzpvertex;
  
  /**
   * Maximum weight for integration
   */
  vector<double> _h0wgt;
  
  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Scalar wavefunction
   */
  mutable ScalarWaveFunction _swave;

  /**
   *  Vector wavefunctions
   */
  mutable vector<VectorWaveFunction> _vwave[2];

private:

  /**
   *  Parameters for the real POWHEG correction
   */
  //@{
  /**
   *  Minimum flavour of quarks to include in the loops
   */
  int _minloop;

  /**
   *  Maximum flavour of quarks to include in the loops
   */
  int _maxloop;

  /**
   *  Option for treatment of the fermion loops
   */
  unsigned int _massopt;
  //@}
};
  
}

#endif /* HERWIG_SMHiggsGGHiggsPPDecayer_H */
