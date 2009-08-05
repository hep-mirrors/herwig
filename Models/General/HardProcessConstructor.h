// -*- C++ -*-
//
// HardProcessConstructor.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HardProcessConstructor_H
#define HERWIG_HardProcessConstructor_H
//
// This is the declaration of the HardProcessConstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "HPDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Utilities/Exception.h"
#include "HardProcessConstructor.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * The HardProcessConstructor is designed to construct the diagrams that are 
 * possible for a given set of incoming and outgoing particles.
 *
 * @see \ref HardProcessConstructorInterfaces "The interfaces"
 * defined for HardProcessConstructor.
 * @see Interfaced
 */
class HardProcessConstructor: public Interfaced {

public:

  /** Set of ParticleData pointers */
  typedef set<tPDPtr> tPDSet;

  /** Vector of HPDiagrams. */
  typedef vector<HPDiagram> HPDVector;

  /** Map of HPDiagrams. */
  typedef multimap<HPDiagram, HPDiagram> HPDMap;

  /** Enumeration for the direction */
  enum direction {incoming, outgoing};

public:

  /**
   * The default constructor.
   */
  HardProcessConstructor();

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

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HardProcessConstructor> initHardProcessConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardProcessConstructor & operator=(const HardProcessConstructor &);

private:

  /**
   * Find if a duplicate pair of particles exists in
   * pair storage
   * @param ppair The ppair to check.
   */
  bool duplicate(tcPDPair ppair) const;

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
  void createTChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex);

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
   * Determine whether the ordering of the outgoing states is the same
   * as the ordering in the matrix elements
   * @param diag The diagram to question
   */
  void fixFSOrder(HPDiagram & diag);

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

  /**
   * Search for a diagram that has already been created
   * @param diagram The diagram to search for
   * @param group The group of diagrams to search through 
   */
  bool duplicate(const HPDiagram & diagram, 
		 const HPDVector & group) const;

  /**
   * Return whether the two ID codes are of the same flavour
   * @param id1 The PDG code of the first particle
   * @param id2 The PDG code of the first particle
   */
  bool sameQuarkFlavour(long id1, long id2) const;
   //@}
  
  /** Functions to set up colour flows and matrix elements. */
  //@{
  /**
   * Assign a diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void assignToCF(HPDiagram & diag);

  /**
   * Assign a $t$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void tChannelCF(HPDiagram & diag);

  /**
   * Assign a $u$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void uChannelCF(HPDiagram & diag);

  /**
   * Assign a $s$-channel diagram to the appropriate colour flow(s).
   * @param diag The diagram to assign
   */
  void sChannelCF(HPDiagram & diag);

  /**
   * Get the correct colour factor matrix.
   * @param extpart Vector of external ParticleData pointers
   * @param ncf Set the number of colourflows.
   */
  vector<DVector> getColourFactors(const tcPDVector & extpart, 
				   unsigned int & ncf) const;
  /**
   * Contruct the classname and object name for the matrix element
   * @param extpart vector containing incoming and outgoing particle data pointers
   * @param objname a string containing the default path of the ME object
   */  
  string MEClassname(const vector<tcPDPtr> & extpart, 
		     string & objname) const;
  //@}
 
private:

  /**
   *  Option for effective vertices
   */
  bool theEffective;

  /**
   * Required initial state particles
   */
  PDVector theIncoming;

  /**
   * Pairs of particles for initial state, ordered by spin or id.
   * If both are of differing spin then lowest is first and if the spin is
   * equal the particle goes first then the anti-particle. This is setup
   * in the doinit() member.
   */
  vector<tcPDPair> theIncPairs;
  
  /**
   * Required final state particles
   */
  PDVector theOutgoing;

  /**
   * Number of incoming particles
   */
  unsigned int theNout;

  /**
   * Pointer to the model being used
   */
  tHwSMPtr theModel;

  /**
   * Number of vertices in the model
   */
  unsigned int theNv;

  /**
   *  The vertices
   */
  vector<VertexBasePtr> theVertices;

  /**
   * Store the configuration of the diagrams
   */
  HPDVector theProcesses;

  /**
   * Whether to include all diagrams or just those with strong
   * coupling in them
   */
  bool theAllDiagrams;

  /**
   * Which types of processes to generate
   */
  unsigned int theProcessOption;

  /**
   * Whether to print the debug information with the matrix 
   * element. This is here solely so it can be passed to 
   * a matrix element that is created here.
   */
  bool theDebug;
  
  /**
   * Pointer to the sub process handler
   */
   tSubHdlPtr theSubProcess;
  
  /**@name The colour factor matrices. */ 
  //@{
  /**
   * The colour factors for the \f$3\,\,\bar{3}\rightarrow 3\,\,\bar{3}\f$.
   */
  vector<DVector> the33bto33b;

  /**
   * The colour factors for the \f$3\,\,\bar{3}^{'}\rightarrow 3\,\,\bar{3}^{'}\f$.
   */
  vector<DVector> the33bpto33bp;

  /**
   * The colour factors for the \f$3\,\,\bar{3}\rightarrow 8\,\,8\f$.
   */
  vector<DVector> the33bto88;

  /**
   * The colour factors for the \f$8\,\,8\rightarrow 8\,\,8\f$.
   */
  vector<DVector> the88to88;
  //@}
};

  /** Exception class indicating setup problem. */
  class HardProcessConstructorError : public Exception {};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HardProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::HardProcessConstructor,1> {
  /** Typedef of the first base class of HardProcessConstructor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HardProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HardProcessConstructor>
  : public ClassTraitsBase<Herwig::HardProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HardProcessConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_HardProcessConstructor_H */
