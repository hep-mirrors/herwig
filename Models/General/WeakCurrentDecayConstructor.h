// -*- C++ -*-
//
// WeakCurrentDecayConstructor.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_WeakCurrentDecayConstructor_H
#define HERWIG_WeakCurrentDecayConstructor_H
//
// This is the declaration of the WeakCurrentDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Decay/General/GeneralCurrentDecayer.fh"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Decay/WeakCurrents/WeakDecayCurrent.h"

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
  inline WeakCurrentDecayConstructor() :
    _theExistingDecayers(0),_init(true),_iteration(5),_points(10000),
    _masscut(5.*GeV) {}

  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const set<PDPtr> & part);

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const { return 2; }

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
  vector<tPDPtr> createModes(const PDPtr inpart,
			    const VertexBasePtr vert,
			    unsigned int ilist,
			    unsigned int iv);

  /**
   * Function to create decayer for specific vertex
   * @param vert Pointer to vertex 
   * @param icol Integer referring to the colmun in _theExistingDecayers
   * @param ivert Integer referring to the row in _theExistingDecayers
   * member variable
   */
  void createDecayer(const VertexBasePtr vert, unsigned int icol,
		     unsigned int ivert);

  /**
   * Create decay mode(s) from given part and decay modes
   * @param inpart pointer to incoming particle
   * @param decays list of allowed interactions
   * @param decayer The decayer responsible for this decay
   */
  void createDecayMode(PDPtr inpart,
		       const tPDVector & decays,
		       map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> decayer);

  /**
   * Set the interfaces on the decayers to initialise them
   * @param name Fullname of the decayer in the EventGenerator
   * including the path
   */
  void initializeDecayers(string name) const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<WeakCurrentDecayConstructor> initWeakCurrentDecayConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeakCurrentDecayConstructor & operator=(const WeakCurrentDecayConstructor &);

private:

  /**
   *  Existing decayers
   */
  vector<vector<map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> > >
  _theExistingDecayers;

  /**
   * Model Pointer
   */
  Ptr<Herwig::StandardModel>::pointer _theModel;

  /**
   * Whether to initialize the decayers or not
   */
  bool _init;
  
  /**
   * Number of iterations if initializing (default 1)
   */
  int _iteration;

  /**
   * Number of points to do in initialization
   */
  int _points;

  /**
   *  Cut-off on the mass difference
   */
  Energy _masscut;


  /**
   *  Particles for the mode
   */
  //@{
  /**
   *  First decay product
   */
  vector<long> _part1;

  /**
   *  Second decay product
   */
  vector<long> _part2;

  /**
   *  Third  decay product
   */
  vector<long> _part3;

  /**
   *  Fourth decay product
   */
  vector<long> _part4;

  /**
   *  Fifth  decay product
   */
  vector<long> _part5;
  //@}

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

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of WeakCurrentDecayConstructor. */
template <>
struct BaseClassTrait<Herwig::WeakCurrentDecayConstructor,1> {
  /** Typedef of the first base class of WeakCurrentDecayConstructor. */
  typedef Herwig::NBodyDecayConstructorBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the WeakCurrentDecayConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::WeakCurrentDecayConstructor>
  : public ClassTraitsBase<Herwig::WeakCurrentDecayConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::WeakCurrentDecayConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_WeakCurrentDecayConstructor_H */
