// -*- C++ -*-
//
// GeneralCurrentDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralCurrentDecayer_H
#define HERWIG_GeneralCurrentDecayer_H
//
// This is the declaration of the GeneralCurrentDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/WeakCurrents/WeakDecayCurrent.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "GeneralCurrentDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

/**
 * Here is the documentation of the GeneralCurrentDecayer class.
 *
 * @see \ref GeneralCurrentDecayerInterfaces "The interfaces"
 * defined for GeneralCurrentDecayer.
 */
class GeneralCurrentDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  GeneralCurrentDecayer() : 
    _maxmass(5.*GeV), _wgtmax(0.) {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
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
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay, MEOption meopt) const = 0;
  
  /**
   * Function to return partial Width
   * @param inpart Pointer to incoming particle data object
   * @param outa Pointer to first outgoing particle data object
   * @param currout Pointer to particles in the current
   */
  virtual Energy partialWidth(tPDPtr inpart, tPDPtr outa,
			      vector<tPDPtr> currout) = 0;
  //@}

  /**
   *  set up the decay
   */
  void setDecayInfo(PDPtr in, PDPtr out, const vector<tPDPtr> & outCurrent,
		    VertexBasePtr vertex, WeakDecayCurrentPtr current,
		    Energy maxmass);

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
   *  The number of the mode
   * @param cc Whether of not this is the charge conjugate of the defined mode
   * @param id The PDG codes of the particles
   */
  int modeNumber(bool & cc, vector<long> id) const;

  /**
   *  Access to the map between the number of the mode and the modes in
   *  the current
   */
  unsigned int mode() const { return _mode; }

  /**
   *  Access to the weak current
   */
  WeakDecayCurrentPtr weakCurrent() const { return _current; }

  /**
   * Get vertex pointer
   * @return a pointer to the vertex
   */
  VertexBasePtr getVertex() const { return _theVertex; }

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<GeneralCurrentDecayer> initGeneralCurrentDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralCurrentDecayer & operator=(const GeneralCurrentDecayer &) = delete;

private:
  
  /**
   * Pointer to vertex set in inheriting class
   */
  VertexBasePtr _theVertex;
  
  /**
   * Incoming particle
   **/
  PDPtr _inpart;

  /**
   * First outgoing particle
   */
  PDPtr _outpart;

  /**
   *  Outgoing particles from the current
   */
  vector<tPDPtr> _currentOut; 

  /**
   * Pointer to the current
   */
  WeakDecayCurrentPtr _current;

  /**
   *  Maximum mass difference
   */
  Energy _maxmass;

  /**
   * mapping of the modes to the currents
   */
  unsigned int _mode;

  /**
   * location of the weights
   */
  int _wgtloc;

  /**
   * the maximum weight
   */
  double _wgtmax;

  /**
   *  The weights for the different channels
   */
  vector<double> _weights;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralCurrentDecayer. */
template <>
struct BaseClassTrait<Herwig::GeneralCurrentDecayer,1> {
  /** Typedef of the first base class of GeneralCurrentDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralCurrentDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralCurrentDecayer>
  : public ClassTraitsBase<Herwig::GeneralCurrentDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GeneralCurrentDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_GeneralCurrentDecayer_H */
