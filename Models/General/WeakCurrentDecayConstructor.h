// -*- C++ -*-
//
// WeakCurrentDecayConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WeakCurrentDecayConstructor_H
#define HERWIG_WeakCurrentDecayConstructor_H
//
// This is the declaration of the WeakCurrentDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Decay/General/GeneralCurrentDecayer.fh"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/WeakCurrents/WeakDecayCurrent.h"
#include "Herwig/Decay/General/GeneralCurrentDecayer.h"
#include "TwoBodyDecay.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the WeakCurrentDecayConstructor class.
 *
 * @see \ref WeakCurrentDecayConstructorInterfaces "The interfaces"
 * defined for WeakCurrentDecayConstructor.
 */
class WeakCurrentDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  WeakCurrentDecayConstructor() : _masscut(5.*GeV) {}
  
  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const set<PDPtr> & part);

  /**
   * Number of outgoing lines. Required for correct ordering (do this one last)
   */
  virtual unsigned int numBodies() const { return 1000; }

  /**
   *  Cut off
   */
  Energy massCut() const { return _masscut;}

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

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

  /** @name Functions to create decayers and decaymodes. */
  //@{
  /**
   * Function to create decays
   * @param inpart Incoming particle 
   * @param vert The vertex to create decays for
   * @param ilist Which list to search
   * @param iv Row number in _theExistingDecayers member
   * @return vector of ParticleData ptrs
   */
  vector<TwoBodyDecay>
  createModes(const PDPtr inpart,const VertexBasePtr vert,
	      unsigned int ilist);

  /**
   * Function to create decayer for specific vertex
   * @param vert Pointer to vertex 
   * @param icol Integer referring to the colmun in _theExistingDecayers
   * @param ivert Integer referring to the row in _theExistingDecayers
   * member variable
   */
  GeneralCurrentDecayerPtr createDecayer(PDPtr in, PDPtr out1,
					 vector<tPDPtr> outCurrent,
					 VertexBasePtr vertex,
					 WeakDecayCurrentPtr current);

  /**
   * Create decay mode(s) from given part and decay modes
   * @param inpart pointer to incoming particle
   * @param decays list of allowed interactions
   * @param decayer The decayer responsible for this decay
   */
  void createDecayMode(vector<TwoBodyDecay> & decays);
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeakCurrentDecayConstructor & operator=(const WeakCurrentDecayConstructor &) = delete;

private:

  /**
   * Model Pointer
   */
  Ptr<Herwig::StandardModel>::pointer _theModel;

  /**
   *  Cut-off on the mass difference
   */
  Energy _masscut;

  /**
   *  Tags for the modes
   */
  vector<string> decayTags_;

  /**
   *  Particles for the mode
   */
  vector<vector<tPDPtr> > particles_;

  /**
   *  Normalisation
   */
  vector<double> _norm;

  /**
   *  The current for the mode
   */
  vector<WeakDecayCurrentPtr> _current;
};

}

#endif /* HERWIG_WeakCurrentDecayConstructor_H */
