// -*- C++ -*-
//
// TwoToTwoProcessConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoToTwoProcessConstructor_H
#define HERWIG_TwoToTwoProcessConstructor_H
//
// This is the declaration of the TwoToTwoProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "ThePEG/Utilities/Exception.h"
#include "TwoToTwoProcessConstructor.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The TwoToTwoProcessConstructor is designed to construct the diagrams that are 
 * possible for a given set of incoming and outgoing particles.
 *
 * @see \ref TwoToTwoProcessConstructorInterfaces "The interfaces"
 * defined for TwoToTwoProcessConstructor.
 * @see HardProcessConstructor
 */
class TwoToTwoProcessConstructor: public HardProcessConstructor {

public:

  /** Set of ParticleData pointers */
  typedef set<tPDPtr> tPDSet;

  /** Map of HPDiagrams. */
  typedef multimap<HPDiagram, HPDiagram> HPDMap;

  /** Enumeration for the direction */
  enum direction {incoming, outgoing};

public:

  /**
   * The default constructor.
   */
  TwoToTwoProcessConstructor();

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

public:

  /**
   * Main function called to start constructing the diagrams for 
   * the 2->2 process
   */
  void constructDiagrams();  

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

  /** @name Standard HardProcessConstructor functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TwoToTwoProcessConstructor> initTwoToTwoProcessConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoToTwoProcessConstructor & operator=(const TwoToTwoProcessConstructor &);

private:

  /** Functions to create the diagrams.*/
  //@{
  /**
   * Given a vertex and 2 particle id's find the possible states
   * that can be the 3rd external particle
   * @param vertex Pointer to the vertex
   * @param part1 id of first particle
   * @param d1 direction of particle one
   * @param part2 id of other particle
   * @param d2 direction of particle two
   * @param d3 required direction of 3rd state (default = outgoing)
   * @return container of third particles
   */
  tPDSet search(VertexBasePtr vertex, long part1, direction d1, 
	       long part2, direction d2, direction d3 = outgoing);    

  /**
   * Given a vertex and 3 particle id's find the possible states
   * that can be the 4th external particle
   * @param vertex Pointer to the vertex
   * @param part1 id of first particle
   * @param d1 direction of particle one
   * @param part2 id of second particle
   * @param d2 direction of particle two
   * @param part3 id of third particle
   * @param d3 direction of particle three
   * @param d4 Required direction of fourth state (default = outgoing)
   * @return container of fourth particles
   */
  tPDSet search(VertexBasePtr vertex, long part1, direction d1, long part2,
	       direction d2, long part3, direction d3, 
	       direction d4 = outgoing);
  
  /**
   * Create the resonance diagrams.
   * @param inpp The incoming pair of particles.
   * @param fs A possible final state.
   * @param vertex The possible interaction vertex for the incoming pair
   */
  void createSChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex);

 /**
   * Create the scattering diagrams.
   * @param inpp The incoming pair of particles.
   * @param fs A possible final state.
   * @param vertex The first vertex
   */
  void createTChannels(tPDPair inpp, long fs, tVertexBasePtr vertex);

  /**
   * Populate the diagram structure
   * @param in Pair of incoming particle id's
   * @param out1 first outgoing particle
   * @param out2 set of second outgoing particles
   * @param inter pointer to particle data for intermediate
   * @param chan the channel type 
   * @param vertices pair of vertices for the diagram
   * @param order The order
   */
  void makeDiagrams(IDPair in, long out1, const tPDSet & out2, PDPtr inter,
		    HPDiagram::Channel chan, VBPair vertices, BPair order);

  /**
   * Create diagrams from 4 point vertices
   * @param parta id of first incoming particle
   * @param partb id of second incoming particle
   * @param partc id of first outgoing particle
   * @param vert pointer to the vertex
   */
  void makeFourPointDiagrams(long parta,long partb,long partc,
			     VertexBasePtr vert);
    
  /**
   * Create the matrix element that will calculate me2() for this
   * process
   * @param process vector of HardPrcoessDiagrams structs that store 
   * the information about the diagrams
   */  
  void createMatrixElement(const HPDVector & process) const;
   //@}
  

  /**
   * Contruct the classname and object name for the matrix element
   * @param extpart vector containing incoming and outgoing particle data pointers
   * @param objname a string containing the default path of the ME object
   */  
  string MEClassname(const vector<tcPDPtr> & extpart, 
		     string & objname) const;
 
private:

  /**
   * Required initial state particles
   */
  PDVector incoming_;

  /**
   * Pairs of particles for initial state, ordered by spin or id.
   * If both are of differing spin then lowest is first and if the spin is
   * equal the particle goes first then the anti-particle. This is setup
   * in the doinit() member.
   */
  vector<tPDPair> incPairs_;
  
  /**
   * Required final state particles
   */
  PDVector outgoing_;

  /**
   * Number of incoming particles
   */
  unsigned int Nout_;

  /**
   * Number of vertices in the model
   */
  unsigned int nv_;

  /**
   *  The vertices
   */
  vector<VertexBasePtr> vertices_;

  /**
   * Store the configuration of the diagrams
   */
  HPDVector processes_;

  /**
   * Whether to include all diagrams or just those with strong
   * coupling in them
   */
  bool allDiagrams_;

  /**
   * Which types of processes to generate
   */
  unsigned int processOption_;

  /**
   *  Option for the scales
   */
  unsigned int scaleChoice_;

  /**
   *  Prefactor for the scale calculation
   */
  double scaleFactor_;

  /**
   *  Option to exclude certain intermediates
   */
  vector<PDPtr> excluded_;

  /**
   *  Option to exclude certain external particles
   */
  vector<PDPtr> excludedExternal_;

  /**
   *  Excluded Vertices
   */
  vector<VertexBasePtr> excludedVertexVector_;

  /**
   *  Excluded Vertices
   */
  set<VertexBasePtr> excludedVertexSet_;
};

  /** Exception class indicating setup problem. */
  class TwoToTwoProcessConstructorError : public Exception {

  public:

    /**
     * Exception for error handling
     * @param str Error message
     * @param sev Severity
     */
    TwoToTwoProcessConstructorError(const string & str, Severity sev) 
      : Exception(str,sev) {}

  };

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TwoToTwoProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::TwoToTwoProcessConstructor,1> {
  /** Typedef of the first base class of TwoToTwoProcessConstructor. */
  typedef Herwig::HardProcessConstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TwoToTwoProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TwoToTwoProcessConstructor>
  : public ClassTraitsBase<Herwig::TwoToTwoProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TwoToTwoProcessConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_TwoToTwoProcessConstructor_H */
