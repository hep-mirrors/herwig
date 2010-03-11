// -*- C++ -*-
//
// GeneralHardME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralHardME_H
#define HERWIG_GeneralHardME_H
//
// This is the declaration of the GeneralHardME class.
//

#include "Herwig++/MatrixElement/HwMEBase.h"
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
 * @see HwMEBase
 */

class GeneralHardME: public HwMEBase {

public:

  /**
   * Convenient typedef for size_type of HPDiagram vector 
   */
  typedef vector<HPDiagram>::size_type HPCount;

public:

  /**
   * The default constructor.
   */
  GeneralHardME();

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
  virtual Energy2 scale() const {
    if(scaleChoice_==0) {
      return sHat();
    }
    else {
      assert( scaleChoice_== 1 );
      Energy2 t = 0.5*(tHat()-meMomenta()[2].mass2());
      Energy2 u = 0.5*(uHat()-meMomenta()[3].mass2());
      Energy2 s = 0.5*sHat();
      return 4.*s*t*u/(s*s+t*t+u*u);
    }
  }

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
   * @param debug Whether to compare the numerical answer to an analytical
   * formula (This is only stored for certain processes. It is intended
   * for quick checks of the matrix elements).
   * @param scaleOption The option of what scale to use
   */
  void setProcessInfo(const vector<HPDiagram> & process,
		      const vector<DVector> & factors,
		      const unsigned int ncf,
		      bool debug, unsigned int scaleOption);

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
   * A debugging function to test the value of me2 against an
   * analytic function. This is to be overidden in an inheriting class.
   * @param x The value of \f$ |\mathcal{M} |^2 \f$
   */
  virtual void debug(double x) const;  

protected:

  /**
   * Access the HPDiagrams that store the required information
   * to create the diagrams
   */
  const vector<HPDiagram> & getProcessInfo() const {
    return theDiagrams;
  }

  /**
   * Return the incoming pair
   * @return Pair of particle ids for the incoming particles
   */
  pair<long, long> getIncoming() const {
    return theIncoming;
  }
  
  /**
   * Return the outgoing pair
   * @return Pair of particle ids for the outgoing particles
   */
  pair<long, long> getOutgoing() const {
    return theOutgoing;
  }
  
  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getColourFactors() const {
    return theColour;
  }

  /**
   * Get the number of diagrams in this process
   */
  HPCount numberOfDiags() const {
    return theNDiags;
  } 
  
  /**
   * Access number of colour flows
   */
  size_t numberOfFlows() const {
    return theNcf;
  }

  /**
   * Whether to print the debug information 
   */
  bool debugME() const {
    return theDebug;
  }

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

  /**
   * Whether to test the value of me2 against the analytical function
   */
  bool theDebug;

  /**
   * The scale chocie
   */
  unsigned int scaleChoice_;

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
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralHardME class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralHardME>
  : public ClassTraitsBase<Herwig::GeneralHardME> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GeneralHardME"; }
};

/** @endcond */

}

#endif /* HERWIG_GeneralHardME_H */
