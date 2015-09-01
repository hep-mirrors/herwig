// -*- C++ -*-
//
// GeneralHardME.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GeneralHardME_H
#define HERWIG_GeneralHardME_H
//
// This is the declaration of the GeneralHardME class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/General/HPDiagram.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"
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

  /**
   *  Enum for the possible colour structures
   */
  enum ColourStructure {UNDEFINED,
                        Colour11to11,Colour11to33bar,Colour11to88,
                        Colour33to33,Colour33barto11,Colour33barto33bar,
                        Colour33barto66bar, Colour33barto6bar6,
			Colour33to61, Colour3bar3barto6bar1,
                        Colour33to16, Colour3bar3barto16bar,
			Colour38to3bar6, Colour38to63bar,
                        Colour33barto18,Colour33barto81,Colour33barto88,
                        Colour38to13,Colour38to31,
                        Colour38to83,Colour38to38,
                        Colour3bar3barto3bar3bar,
                        Colour3bar8to13bar,Colour3bar8to3bar1,
                        Colour3bar8to83bar,Colour3bar8to3bar8,
                        Colour88to11,Colour88to33bar,
                        Colour88to66bar,Colour88to88,
                        Colour88to18,Colour88to81};

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
      return scaleFactor_*sHat();
    }
    else if(scaleChoice_==1) {
      Energy2 mbar = 0.5*(meMomenta()[2].mass2()+meMomenta()[3].mass2());
      Energy2 t = 0.5*(tHat()-mbar);
      Energy2 u = 0.5*(uHat()-mbar);
      Energy2 s = 0.5*sHat();
      return scaleFactor_*4.*s*t*u/(s*s+t*t+u*u);
    }
    else if(scaleChoice_ ==2) {
      Energy2 scale1 = meMomenta()[2].mass2()+meMomenta()[2].perp2();
      Energy2 scale2 = meMomenta()[3].mass2()+meMomenta()[3].perp2();
      return scaleFactor_*max(scale1,scale2);
    }
    else assert(false);
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
  colourGeometries(tcDiagPtr diag) const;
  //@}

  /**
   * Set the diagrams and matrix of colour factors. 
   * @param process vector of MEDiagram with information that 
   * will allow the diagrams to be created in the specific matrix element
   * @param colour The colour structure for the process
   * @param debug Whether to compare the numerical answer to an analytical
   * formula (This is only stored for certain processes. It is intended
   * for quick checks of the matrix elements).
   * @param scaleOption The option of what scale to use
   * @param scaleFactor The prefactor for the scale
   */
  void setProcessInfo(const vector<HPDiagram> & process,
		      ColourStructure colour, bool debug, 
		      unsigned int scaleOption,
		      double scaleFactor);

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   * A debugging function to test the value of me2 against an
   * analytic function. This is to be overidden in an inheriting class.
   */
  virtual void debug(double ) const {}  

protected:

  /**
   * Access the HPDiagrams that store the required information
   * to create the diagrams
   */
  const vector<HPDiagram> & getProcessInfo() const {
    return diagrams_;
  }

  /**
   * Return the incoming pair
   * @return Pair of particle ids for the incoming particles
   */
  pair<long, long> getIncoming() const {
    return incoming_;
  }
  
  /**
   * Return the outgoing pair
   * @return Pair of particle ids for the outgoing particles
   */
  pair<long, long> getOutgoing() const {
    return outgoing_;
  }
  
  /**
   * Return the matrix of colour factors 
   */
  const vector<DVector> & getColourFactors() const {
    return colour_;
  }

  /**
   * Get the number of diagrams in this process
   */
  HPCount numberOfDiags() const {
    return numberOfDiagrams_;
  } 
  
  /**
   * Access number of colour flows
   */
  size_t numberOfFlows() const {
    return numberOfFlows_;
  }

  /**
   * Whether to print the debug information 
   */
  bool debugME() const {
    return debug_;
  }

  /**
   *  Set/Get Info on the selected diagram and colour flow
   */
  //@{
  /**
   * Colour flow
   */
  unsigned int colourFlow() const {return flow_;}

  /**
   * Colour flow
   */
  void colourFlow(unsigned int flow) const {flow_=flow;}

  /**
   * Diagram
   */
  unsigned int diagram() const {return diagram_;}

  /**
   * Diagram
   */
  void diagram(unsigned int diag) const {diagram_=diag;}
  //@}

  /**
   *  Calculate weight and select colour flow
   */
  double selectColourFlow(vector<double> & flow,
			  vector<double> & me,double average) const;

  /**
   *  Access to the colour flow matrix element
   */
  vector<ProductionMatrixElement> & flowME() const {
    return flowME_;
  }

  /**
   *  Access to the diagram matrix element
   */
  vector<ProductionMatrixElement> & diagramME() const {
    return diagramME_;
  }

  /**
   *  Access to the colour structure
   */
  ColourStructure colour() const {return colourStructure_;}

  /**
   *  Extract the paricles from the subprocess
   */
  ParticleVector hardParticles(tSubProPtr subp) {
    ParticleVector output(4);
    output[0] = subp->incoming().first; 
    output[1] = subp->incoming().second;
    output[2] = subp->outgoing()[0]; 
    output[3] = subp->outgoing()[1];    
    //ensure particle ordering is the same as it was when
    //the diagrams were created
    if( output[0]->id() != getIncoming().first )
      swap(output[0], output[1]);
    if( output[2]->id() != getOutgoing().first )
      swap(output[2], output[3]);
    // return answer
    return output;
  }

  /**
   *  Set the rescaled momenta
   */
  void setRescaledMomenta(const ParticleVector & external) {
    cPDVector data(4);
    vector<Lorentz5Momentum> momenta(4);
    for( size_t i = 0; i < 4; ++i ) {
      data[i] = external[i]->dataPtr();
      momenta[i] = external[i]->momentum();
    }
    rescaleMomenta(momenta, data);
  }

  /**
   *  Create the vertes
   */
  void createVertex(ProductionMatrixElement & me,
		    ParticleVector & external) {
    HardVertexPtr hardvertex = new_ptr(HardVertex());
    hardvertex->ME(me);
    for(ParticleVector::size_type i = 0; i < 4; ++i) {
      tSpinPtr spin = external[i]->spinInfo();
      if(i<2) {
	tcPolarizedBeamPDPtr beam = 
	  dynamic_ptr_cast<tcPolarizedBeamPDPtr>(external[i]->dataPtr());
	if(beam) spin->rhoMatrix() = beam->rhoMatrix();
      }
      spin->productionVertex(hardvertex);
    }
  }

  /**
   *  Initialize the storage of the helicity matrix elements
   */
  void initializeMatrixElements(PDT::Spin  in1, PDT::Spin in2,
				PDT::Spin out1, PDT::Spin out2) {
    flowME().resize(numberOfFlows(),
		    ProductionMatrixElement(in1,in2,out1,out2));
    diagramME().resize(numberOfDiags(),
		       ProductionMatrixElement(in1,in2,out1,out2));
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
   *  External particles
   */
  //@{
  /**
   * Store incoming particles
   */
  pair<long, long> incoming_;
  
  /**
   * Store the outgoing particles
   */
  pair<long, long> outgoing_;
  //@}

  /**
   *  Diagrams
   */
  //@{
  /**
   * Store all diagrams as a vector of structures
   */
  vector<HPDiagram> diagrams_;

  /**
   * Store the number of diagrams for fast retrieval
   */
  HPCount numberOfDiagrams_;
  //@}

  /**
   *  Colour information
   */
  //@{
  /**
   *  The colour structure
   */
  ColourStructure colourStructure_;

  /**
   * Store colour factors for ME calc.
   */
  vector<DVector> colour_;

  /**
   * The number of colourflows.
   */
  unsigned int numberOfFlows_;
  //@}

  /**
   * Whether to test the value of me2 against the analytical function
   */
  bool debug_;

  /**
   * The scale chocie
   */
  unsigned int scaleChoice_;

  /**
   *  The scale factor
   */
  double scaleFactor_;

  /**
   *  Info on the selected diagram and colour flow
   */
  //@{
  /**
   * Colour flow
   */
  mutable unsigned int flow_;

  /**
   * Diagram
   */
  mutable unsigned int diagram_;
  //@}

  /**
   *  Storage of the matrix elements
   */
  //@{
  /**
   *  Matrix elements for the different colour flows
   */
  mutable vector<ProductionMatrixElement> flowME_;

  /**
   *  Matrix elements for the different Feynman diagrams
   */
  mutable vector<ProductionMatrixElement> diagramME_;
  //@}

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
