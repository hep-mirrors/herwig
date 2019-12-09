// -*- C++ -*-
//
// TwoBodyDecayConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/ShowerAlpha.h"
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
  TwoBodyDecayConstructor() : inter_(ShowerInteraction::Both) {
    radiationVertices_[ShowerInteraction::QCD] = map<tPDPtr,VertexBasePtr>();
    radiationVertices_[ShowerInteraction::QED] = map<tPDPtr,VertexBasePtr>();
  }

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoBodyDecayConstructor & operator=(const TwoBodyDecayConstructor &) = delete;

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
  void createModes(tPDPtr inpart, VertexBasePtr vert,
		   unsigned int ilist,
		   multiset<TwoBodyDecay> & modes);

  /**
   * Function to create decayer for specific vertex
   * @param decay decay mode for this decay
   * member variable
   */
  GeneralTwoBodyDecayerPtr createDecayer(TwoBodyDecay decay,
					 vector<VertexBasePtr> );

  /**
   * Create decay mode(s) from given part and decay modes
   * @param decays The vector of decay modes
   * @param decayer The decayer responsible for this decay
   */
  void createDecayMode(multiset<TwoBodyDecay> & decays);
  //@}

  /**
   * Get the vertex for QED/QCD radiation
   */
  VertexBasePtr radiationVertex(tPDPtr particle,ShowerInteraction inter,
				tPDPair children = tPDPair ());

private:

  /**
   *  Map of particles and the vertices which generate their QCD
   *  radiation
   */
  map<ShowerInteraction,map<tPDPtr,VertexBasePtr> > radiationVertices_;

  /**
   *  Default choice for the strong coupling object for hard QCD radiation
   */
  ShowerAlphaPtr  alphaQCD_;

  /**
   *  Default choice for the strong coupling object for hard QED radiation
   */
  ShowerAlphaPtr  alphaQED_;

  /**
   *  Which type of corrections to the decays to include
   */
  ShowerInteraction inter_; 
  
};
  
}

#endif /* HERWIG_TwoBodyDecayConstructor_H */
