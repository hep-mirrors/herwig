// -*- C++ -*-
#ifndef HERWIG_ResonantProcessConstructor_H
#define HERWIG_ResonantProcessConstructor_H
//
// This is the declaration of the ResonantProcessConstructor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "HPDiagram.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/SubProcessHandler.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
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
 * @see Interfaced
 */
class ResonantProcessConstructor: public Interfaced {

public:

  /** Set of ParticleData pointers */
  typedef set<PDPtr> PDSet;

  /** Vector of HPDiagrams. */
  typedef vector<HPDiagram> HPDVector;

  /** Nested vector of doubles. */
  typedef vector<vector<double> > CFMatrix;

  /** Enumeration for the direction */
  enum direction {incoming, outgoing};

public:

  /**
   * The default constructor.
   */
  inline ResonantProcessConstructor();

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
  void constructResonances() ;

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
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
  void makeResonantDiagrams(IDPair in, PDPtr offshell, long outa, 
			    const PDSet & out, VBPair vertices);
  
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
  PDSet search(VertexBasePtr vertex, long part1, direction d1, 
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
  string MEClassname(const tcPDVector & extpart, string & objname) const;

  /**
   * Return colour factor for given process
   */
  CFMatrix colourFactor(const tcPDVector & extpart) const;

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
   * Storage for the required intermediate particles
   */
  vector<PDPtr> theIncoming;

  /**
   * Storage for the required intermediate particles
   */
  vector<PDPtr> theIntermediates;

  /**
   * Storage for the required intermediate particles
   */
  vector<PDPtr> theOutgoing;  
  
  /**
   * A pointer to the model being used
   */
  HwSMPtr theModel;

  /**
   * Storage for the diagrams
   */
  vector<HPDiagram> theDiagrams;

  /**
   * Store a pointer to the subprocess handler
   */
  SubHdlPtr theSubProcess; 
};

  /** Exception class indicating setup problem. */
  class RPConstructorError : public Exception {};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/// \if TRAITSPECIALIZATIONS

/** This template specialization informs ThePEG about the
 *  base classes of ResonantProcessConstructor. */
template <>
struct BaseClassTrait<Herwig::ResonantProcessConstructor,1> {
  /** Typedef of the first base class of ResonantProcessConstructor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ResonantProcessConstructor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ResonantProcessConstructor>
  : public ClassTraitsBase<Herwig::ResonantProcessConstructor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::ResonantProcessConstructor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ResonantProcessConstructor is implemented. It may also include several, space-separated,
   * libraries if the class ResonantProcessConstructor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwModelGenerator.so"; }
};

/// \endif

}

#include "ResonantProcessConstructor.icc"

#endif /* HERWIG_ResonantProcessConstructor_H */
