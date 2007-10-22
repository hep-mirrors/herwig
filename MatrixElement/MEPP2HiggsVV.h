// -*- C++ -*-
#ifndef HERWIG_MEPP2HiggsVV_H
#define HERWIG_MEPP2HiggsVV_H
//
// This is the declaration of the MEPP2HiggsVV class.
//
#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Models/General/SVVLoopVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ProductionMatrixElement.h"
#include "MEPP2HiggsVV.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2HiggsVV class implements the matrix element for the process 
 * pp->Higgs->VV (V = W,Z) with the corrected Higgs width treatment 
 * (see details in hep-ph/9505211).
 *
 * @see \ref MEPP2HiggsVVInterfaces "The interfaces"
 * defined for MEPP2HiggsVV.
 */
class MEPP2HiggsVV: public ME2to2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MEPP2HiggsVV();

  /**
   * The copy constructor.
   */
  inline MEPP2HiggsVV(const MEPP2HiggsVV &);

  /**
   * The destructor.
   */
  virtual ~MEPP2HiggsVV();
  //@}

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
//  virtual CrossSection dSigHatDR() const;

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
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

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
  inline virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual bool generateKinematics(const double * r);
  //@}

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
//  inline virtual void rebind(const TranslationMap & trans)
//    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
//  inline virtual IVector getReferences();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2HiggsVV> initMEPP2HiggsVV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsVV & operator=(const MEPP2HiggsVV &);

  /**
   *  Members to return the matrix elements for the different subprocesses
   */
  //@{

  /**
   * Calculates the matrix element for the process g,g->h (via quark loops)
   * @param g1 a vector of wave functions of the first incoming gluon
   * @param g2 a vector of wave functions of the second incoming gluon
   * @param v1 a vector of the first outgoing boson wave function
   * @param v2 a vector of the second outgoing boson wave function
   * @param calc Whether or not to calculate the matrix element for spin correlations
   * @return the amlitude value.
   */
  double ggME(vector<VectorWaveFunction> g1, 
              vector<VectorWaveFunction> g2, 
              vector<VectorWaveFunction> v1, 
              vector<VectorWaveFunction> v2, 
              bool calc) const;

  /**
   * Calculates the matrix element for the process q,qbar->h
   * @param fin a vector of incoming quark spinors
   * @param ain a vector of incoming anti-quark spinors
   * @param v1 a vector of the first outgoing boson wave function
   * @param v2 a vector of the second outgoing boson wave function
   * @param calc Whether or not to calculate the matrix element for spin correlations
   * @return the amlitude value.
   */
  double qqME(vector<SpinorWaveFunction>    fin, 
              vector<SpinorBarWaveFunction> ain, 
              vector<VectorWaveFunction> v1, 
              vector<VectorWaveFunction> v2, 
              bool calc) const;

private:

  /**
   *  Use the mass generator for the line shape
   */
  unsigned int shapeopt;

  /**
   *  The processes to be included
   */
  unsigned int processopt;

  /**
   *  Minimum flavour of incoming quarks
   */
  unsigned int minflavouropt;

  /**
   *  Maximum flavour of incoming quarks
   */
  unsigned int maxflavouropt;

  /**
   *  Maximum flavour of incoming quarks
   */
  unsigned int bosonopt;

  /**
   *  Storage of the diagram weights for the \f$gg\to Hg\f$ subprocess
   */
  mutable double diagwgt[3];

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   * Pointer to the H->2gluons vertex
   */
  SVVLoopVertexPtr hggvertex;

  /**
   * Pointer to the H->2gamma vertex
   */
  SVVLoopVertexPtr hppvertex;

  /**
   * Pointer to the H->Z,Z/W-,W+ vertex
   */
  VVSVertexPtr vvhvertex;

  /**
   * Pointer to the H->f,f vertex
   */
  FFSVertexPtr ffhvertex;

  /**
   * Pointer to the Standard Model instance used in the class
   */
  tcHwSMPtr theSM;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2HiggsVV. */
template <>
struct BaseClassTrait<Herwig::MEPP2HiggsVV,1> {
  /** Typedef of the first base class of MEPP2HiggsVV. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2HiggsVV class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2HiggsVV>
  : public ClassTraitsBase<Herwig::MEPP2HiggsVV> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2HiggsVV"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2HiggsVV is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2HiggsVV depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwME.so"; }
};

/** @endcond */

}

#include "MEPP2HiggsVV.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEPP2HiggsVV.tcc"
#endif

#endif /* HERWIG_MEPP2HiggsVV_H */
