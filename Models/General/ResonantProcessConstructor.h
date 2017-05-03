// -*- C++ -*-
//
// ResonantProcessConstructor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ResonantProcessConstructor_H
#define HERWIG_ResonantProcessConstructor_H
//
// This is the declaration of the ResonantProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "ThePEG/Utilities/Exception.h"
#include "ResonantProcessConstructor.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * This class is designed to construct the diagrams for resonant processes
 * using a provdided set of particles as interemdiates.
 *
 * @see \ref ResonantProcessConstructorInterfaces "The interfaces"
 * defined for ResonantProcessConstructor.
 * @see HardProcessConstructor
 */
class ResonantProcessConstructor: public HardProcessConstructor {

public:

  /** Set of ParticleData pointers */
  typedef set<tPDPtr> tPDSet;

  /** Nested vector of doubles. */
  typedef vector<vector<double> > CFMatrix;

  /** Enumeration for the direction */
  enum direction {incoming, outgoing};

public:

  /**
   * The default constructor.
   */
  ResonantProcessConstructor() :
    processOption_(0), scaleChoice_(1), scaleFactor_(1.),
    incoming_(0), intermediates_(0),
    outgoing_(0), diagrams_(0)
  {}

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

  /**
   * The main function to create the resonant diagrams
   */
  void constructDiagrams() ;

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
   * Utility function to help second vertex
   */
  void constructVertex2(IDPair in, VertexBasePtr vertex, 
			PDPtr partc);
  
  /**
   * Function to create the appropriate diagrams
   */
  void makeResonantDiagram(IDPair in, PDPtr offshell, long outa, 
			   long outb, VBPair vertices);
  
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
   * Return the pair of outgoing particles from the list
   */
  IDPair find(long part, const PDVector & out) const;

  /**
   * Create a matrix element from the given resonant process diagram
   */
  void createMatrixElement(const HPDiagram & diag) const;

  /**
   * Create the correct classname and objectname for a matrix element 
   */
  string MEClassname(const tcPDVector & extpart, tcPDPtr inter,
		     string & objname) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ResonantProcessConstructor> 
  initResonantProcessConstructor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ResonantProcessConstructor & operator=(const ResonantProcessConstructor &);

private:

  /**
   * Which types of processes to generate
   */
  unsigned int processOption_;

  /**
   *  Scale choice
   */
  unsigned int scaleChoice_;

  /**
   *  Prefactor for the scale calculation
   */
  double scaleFactor_;
  
  /**
   * Storage for the required intermediate particles
   */
  vector<PDPtr> incoming_;

  /**
   * Storage for the required intermediate particles
   */
  vector<PDPtr> intermediates_;

  /**
   * Storage for the required intermediate particles
   */
  vector<PDPtr> outgoing_; 

  /**
   * Storage for the diagrams
   */
  vector<HPDiagram> diagrams_;

};

  /** Exception class indicating setup problem. */
  class RPConstructorError : public Exception {};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ResonantProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::ResonantProcessConstructor,1> {
  /** Typedef of the first base class of ResonantProcessConstructor. */
  typedef Herwig::HardProcessConstructor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ResonantProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ResonantProcessConstructor>
  : public ClassTraitsBase<Herwig::ResonantProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ResonantProcessConstructor"; }
};

/** @endcond */

}

#endif /* HERWIG_ResonantProcessConstructor_H */
