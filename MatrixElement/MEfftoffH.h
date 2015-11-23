// -*- C++ -*-
#ifndef HERWIG_MEfftoffH_H
#define HERWIG_MEfftoffH_H
//
// This is the declaration of the MEfftoffH class.
//

#include "HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEfftoffH class is the base class for vector boson fusion type
 * processes in Herwig.
 *
 * @see \ref MEfftoffHInterfaces "The interfaces"
 * defined for MEfftoffH.
 */
class MEfftoffH: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEfftoffH() : _shapeopt(2), _maxflavour(5), _minflavour(1), _process(0), 
		_mh(), _wh(), _swap(false) {}
  
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
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
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

  /**
   * Matrix element for \f$ff\to h^0\to ff h^0\f$.
   * @param f1  Spinors for first  incoming fermion
   * @param f2  Spinors for second incoming fermion
   * @param a1  Spinors for first  outgoing fermion
   * @param a2  Spinors for second outgoing fermion
   * @param swap1 Whether or not to swap the order for the first fermion line
   * @param swap2 Whether or not to swap the order for the second fermion line
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction> & f1 ,
		    vector<SpinorWaveFunction> & f2 ,
		    vector<SpinorBarWaveFunction> & a1,
		    vector<SpinorBarWaveFunction> & a2,
		    bool swap1, bool swap2,
		    bool me) const;

  /**
   *  Access to the vector ParticleData objects
   */
  //@{
  /**
   *  Access to the \f$W^+\f$ data
   */ 
  PDPtr WPlus() const {return _wplus;}

  /**
   *  Access to the \f$W^-\f$ data
   */ 
  PDPtr WMinus() const {return _wminus;}

  /**
   *  Access to the \f$Z^0\f$ data
   */ 
  PDPtr Z0() const {return _z0;}

  /**
   *  Access to the Higgs boson
   */
  PDPtr higgs() const {return _higgs;}

  /**
   *  Set the Higgs boson
   */
  void higgs(PDPtr in) {_higgs=in;}
  //@}

  /**
   *  Set the pointer to the vector-vector-Higgs vertex
   */
  void setWWHVertex(AbstractVVSVertexPtr in) {
    _vertexWWH = in;
  }

  /**
   *  Set the line shape treatment
   */
  void lineShape(unsigned int in) {_shapeopt=in;}

  /**
   *  Which process to generate
   */
  unsigned int process() const {return _process;}

  /**
   *  Whether momenta are swapped
   */
  bool swapOrder() {return _swap;}

  /**
   *  Which process to generate
   */
  void process(unsigned int in) {_process = in;}

  /**
   *  Maximum flavour of the incoming partons
   */
  unsigned int maxFlavour() const {return _maxflavour;}

  /**
   *  Minimum flavour of the incoming partons
   */
  unsigned int minFlavour() const {return _minflavour;}

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
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<MEfftoffH> initMEfftoffH;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfftoffH & operator=(const MEfftoffH &);

private:

  /**
   * Defines the Higgs resonance shape
   */
  unsigned int _shapeopt;

  /**
   *  Maximum flavour of the quarks involved 
   */
  unsigned int _maxflavour;

  /**
   *  Minimum flavour of the quarks involved 
   */
  unsigned int _minflavour;

  /**
   *  Whether to include $WW$ and $ZZ$ processes or both
   */
  unsigned int _process;

  /**
   *  The intermediate vector bosons
   */
  //@{
  /**
   *  \f$W^+\f$
   */
  PDPtr _wplus;

  /**
   *  \f$W^-\f$
   */
  PDPtr _wminus;

  /**
   *  \f$Z^0\f$
   */
  PDPtr _z0;

  /**
   *  Higgs boson
   */
  PDPtr _higgs;
  //@}

  /**
   *  The vertices for the calculation of the matrix element
   */
  //@{
  /**
   *  Vertex for fermion-fermion-W
   */
  AbstractFFVVertexPtr _vertexFFW;

  /**
   *  Vertex for fermion-fermion-Z
   */
  AbstractFFVVertexPtr _vertexFFZ;

  /**
   *  Vertex for vector-vector-Higgs
   */
  AbstractVVSVertexPtr _vertexWWH;
  //@}

  /**
   *  On-shell mass for the higgs
   */
  Energy _mh;

  /**
   *  On-shell width for the higgs
   */
  Energy _wh;

  /**
   *  The mass generator for the Higgs
   */
  GenericMassGeneratorPtr _hmass;


  /**
   * Matrix element for spin correlations
   */
  mutable ProductionMatrixElement _me;

  /**
   *  if order swaped
   */
  bool _swap;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEfftoffH. */
template <>
struct BaseClassTrait<Herwig::MEfftoffH,1> {
  /** Typedef of the first base class of MEfftoffH. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEfftoffH class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEfftoffH>
  : public ClassTraitsBase<Herwig::MEfftoffH> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEfftoffH"; }
};

/** @endcond */

}

#endif /* HERWIG_MEfftoffH_H */
