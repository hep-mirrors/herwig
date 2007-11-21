// -*- C++ -*-
//
// SMTopDecayer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMTopDecayer_H
#define HERWIG_SMTopDecayer_H
//
// This is the declaration of the SMTopDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "SMTopDecayer.fh"

namespace Herwig {
  using namespace ThePEG;
  using namespace ThePEG::Helicity;
  
/**
 * \ingroup Decay
 *
 * The SMTopDecayer performs decays of the top quark into
 * the bottom quark and qqbar pairs or to the bottom quark and lepton 
 * neutrino pairs via W boson exchange.
 */
class SMTopDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  SMTopDecayer();

  /**
   * Which of the possible decays is required
   */
  virtual int modeNumber(bool & , tcPDPtr , const PDVector & ) const {return -1;}

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
   * Initialize this object after the setup phase before saving and
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
  static ClassDescription<SMTopDecayer> initSMTopDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMTopDecayer & operator=(const SMTopDecayer &);
  
  /**
   *Pointer to the W vertex
   */
  FFVVertexPtr _wvertex;
  
  /**
   * Max weight for integration
   */
  //@{   
  /**
   * Weight \f$W\to q\bar{q}'\f$
   */
  vector<double> _wquarkwgt;
  
  /**
   * Weight \f$W\to \ell \nu\f$
   */
  vector<double> _wleptonwgt;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMTopDecayer. */
template <>
struct BaseClassTrait<Herwig::SMTopDecayer,1> {
  /** Typedef of the first base class of SMTopDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMTopDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMTopDecayer>
  : public ClassTraitsBase<Herwig::SMTopDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMTopDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SMTopDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwPerturbativeDecay.so"; }
};

/** @endcond */

}

#include "SMTopDecayer.icc"

#endif /* HERWIG_SMTopDecayer_H */
