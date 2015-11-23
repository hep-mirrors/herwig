// -*- C++ -*-
#ifndef HERWIG_MEfftoVH_H
#define HERWIG_MEfftoVH_H
//
// This is the declaration of the MEfftoVH class.
//

#include "DrellYanBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The MEfftoVH class is the base class for \f$f\bar{f}\to VH\f$ processes. 
 * This base class handles the phase-space integration while
 * the inheriting classes implement the matrix element
 *
 * @see \ref MEfftoVHInterfaces "The interfaces"
 * defined for MEfftoVH.
 */
class MEfftoVH: public DrellYanBase {

public:

  /**
   * The default constructor.
   */
  MEfftoVH() : _shapeopt(2), _maxflavour(5), _mh(), _wh() {}

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
   * Matrix element for \f$f\bar{f}\to V h^0\to f'\bar{f'} h^0\f$.
   * @param fin  Spinors for incoming fermion
   * @param ain  Spinors for incoming antifermion
   * @param fout Spinors for incoming fermion
   * @param aout Spinors for incoming antifermion
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction> & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<SpinorBarWaveFunction> & fout,
		 vector<SpinorWaveFunction>    & aout,
		 bool me) const;

  /**
   *  Access to the vector ParticleData objects
   */
  //@{
  /**
   *  Access to the \f$W^+\f$ data
   */ 
  PDPtr WPlus() const { return _wplus; }

  /**
   *  Access to the \f$W^-\f$ data
   */ 
  PDPtr WMinus() const { return _wminus; }

  /**
   *  Access to the \f$Z^0\f$ data
   */ 
  PDPtr Z0() const { return _z0; }

  /**
   *  Access to the higgs data
   */ 
  PDPtr higgs() const { return _higgs; }

  /**
   *  Set the higgs data
   */ 
  void higgs(PDPtr in) {_higgs =in;}
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
   *  Maximum flavour of the incoming partons
   */
  unsigned int maxFlavour() const {return _maxflavour;}

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
  static AbstractClassDescription<MEfftoVH> initMEfftoVH;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfftoVH & operator=(const MEfftoVH &);

private:

  /**
   * Defines the Higgs resonance shape
   */
  unsigned int _shapeopt;

  /**
   *  The allowed flavours of the incoming quarks
   */
  unsigned int _maxflavour;

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
   *  The higgs bosom
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
  ProductionMatrixElement _me;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEfftoVH. */
template <>
struct BaseClassTrait<Herwig::MEfftoVH,1> {
  /** Typedef of the first base class of MEfftoVH. */
  typedef Herwig::DrellYanBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEfftoVH class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEfftoVH>
  : public ClassTraitsBase<Herwig::MEfftoVH> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEfftoVH"; }
};

/** @endcond */

}

#endif /* HERWIG_MEfftoVH_H */
