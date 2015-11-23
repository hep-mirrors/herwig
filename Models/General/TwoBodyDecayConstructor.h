// -*- C++ -*-
//
// TwoBodyDecayConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoBodyDecayConstructor_H
#define HERWIG_TwoBodyDecayConstructor_H
//
// This is the declaration of the TwoBodyDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.fh"
#include "TwoBodyDecay.h"

namespace Herwig {
using namespace ThePEG;

using Helicity::VertexBasePtr;
using Helicity::tVertexBasePtr;

/**
 * The TwoBodyDecayConstructor class inherits from the dummy base class
 * NBodyDecayConstructorBase and implements the necessary functions in
 * order to create the 2 body decay modes for a given set of vertices
 * stored in a Model class.
 *
 * @see \ref TwoBodyDecayConstructorInterfaces "The interfaces"
 * defined for TwoBodyDecayConstructor.
 * @see NBodyDecayConstructor
 **/
class TwoBodyDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  TwoBodyDecayConstructor() : showerAlpha_("/Herwig/Shower/AlphaQCD") {}

  /**
   * Function used to determine allowed decaymodes
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const set<PDPtr> & part);

  /**
   * Number of outgoing lines. Required for correct ordering.
   */
  virtual unsigned int numBodies() const { return 2; }

  
public:

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<TwoBodyDecayConstructor> initTwoBodyDecayConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoBodyDecayConstructor & operator=(const TwoBodyDecayConstructor &);

private:
  
  /** @name Functions to create decayers and decaymodes. */
  //@{
  /**
   * Function to create decays
   * @param inpart Incoming particle 
   * @param vert The vertex to create decays for
   * @param ilist Which list to search
   * @param iv Row number in _theExistingDecayers member
   * @return A vector a decay modes
   */
  set<TwoBodyDecay> createModes(tPDPtr inpart, VertexBasePtr vert,
				unsigned int ilist);

  /**
   * Function to create decayer for specific vertex
   * @param decay decay mode for this decay
   * member variable
   */
  GeneralTwoBodyDecayerPtr createDecayer(TwoBodyDecay decay);

  /**
   * Create decay mode(s) from given part and decay modes
   * @param decays The vector of decay modes
   * @param decayer The decayer responsible for this decay
   */
  void createDecayMode(set<TwoBodyDecay> & decays);
  //@}

  /**
   * Get the vertex for QCD radiation
   */
  VertexBasePtr radiationVertex(tPDPtr particle,tPDPair children = tPDPair ());

private:

  /**
   *  Map of particles and the vertices which generate their QCD
   *  radiation
   */
  map<tPDPtr,VertexBasePtr> radiationVertices_;

  /**
   *  Default choice for the strong coupling object for hard radiation
   */
  string showerAlpha_;
};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TwoBodyDecayConstructor. */
template <>
struct BaseClassTrait<Herwig::TwoBodyDecayConstructor,1> {
  /** Typedef of the first base class of TwoBodyDecayConstructor. */
  typedef Herwig::NBodyDecayConstructorBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TwoBodyDecayConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TwoBodyDecayConstructor>
  : public ClassTraitsBase<Herwig::TwoBodyDecayConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TwoBodyDecayConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_TwoBodyDecayConstructor_H */
