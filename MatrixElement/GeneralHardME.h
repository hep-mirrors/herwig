
// -*- C++ -*-
#ifndef HERWIG_GeneralHardME_H
#define HERWIG_GeneralHardME_H
//
// This is the declaration of the GeneralHardME class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/General/HPDiagram.h"
#include "GeneralHardME.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::VertexBasePtr;

/**
 * This defines the GeneralHardME class that is designed to serve as a 
 * base class for matrix elements of specific spin structures when those
 * structures are created by a a general model, i.e. a SUSY production 
 * ME. It stores a vector of diagram structures that contain the required
 * to calculate the matrix element.
 *
 * @see ME2to2Base
 */

class GeneralHardME: public ME2to2Base {

public:

  /**
   * Convenient typedef for size_type of HPDiagram vector 
   */
  typedef vector<HPDiagram>::size_type HPCount;

public:

  /**
   * The default constructor.
   */
  inline GeneralHardME();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const = 0;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual inline Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> 
  diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const = 0;
  //@}

  /**
   * Set the diagrams and matrix of colour factors. 
   * @param process vector of MEDiagram with information that 
   * will allow the diagrams to be created in the specific matrix element
   * @param factors
   * @param ncf Number of colour flows
   */
  inline void setProcessInfo(const vector<HPDiagram> & process,
			     const vector<DVector> & factors,
			     const unsigned int ncf);

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

  /**
   * Access the HPDiagrams that store the required information
   * to create the diagrams
   */
  inline const vector<HPDiagram> & getProcessInfo() const; 

  /**
   * Return the incoming pair
   * @return Pair of particle ids for the incoming particles
   */
  inline pair<long, long> getIncoming() const;

  /**
   * Return the outgoing pair
   * @return Pair of particle ids for the outgoing particles
   */
  inline pair<long, long> getOutgoing() const;
  
  /**
   * Return the matrix of colour factors 
   */
  inline const vector<DVector> & getColourFactors() const;

  /**
   * Get the number of diagrams in this process
   */
  inline const HPCount numberOfDiags() const;
  
  /**
   * Access number of colour flows
   */
  inline const size_t numberOfFlows() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<GeneralHardME> initGeneralHardME;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralHardME & operator=(const GeneralHardME &);

private:
  
  /**
   * Store incoming particles
   */
  pair<long, long> theIncoming;
  
  /**
   * Store the outgoing particles
   */
  pair<long, long> theOutgoing;

  /**
   * Store all diagrams as a vector of structures
   */
  vector<HPDiagram> theDiagrams;

  /**
   * Store the number of diagrams for fast retrieval
   */
  HPCount theNDiags;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> theColour;

  /**
   * The number of colourflows.
   */
  unsigned int theNcf;
  
};

  /** Exception class to indicate a problem has occurred with setting
   up to matrix element.*/
  class MEException : public Exception {};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralHardME. */
template <>
struct BaseClassTrait<Herwig::GeneralHardME,1> {
  /** Typedef of the first base class of GeneralHardME. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralHardME class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralHardME>
  : public ClassTraitsBase<Herwig::GeneralHardME> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GeneralHardME"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GeneralHardME is implemented. It may also include several, space-separated,
   * libraries if the class GeneralHardME depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwGeneralME.so"; }
};

/** @endcond */

}

#include "GeneralHardME.icc"

#endif /* HERWIG_GeneralHardME_H */
