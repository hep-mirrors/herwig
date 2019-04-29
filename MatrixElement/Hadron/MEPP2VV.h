// -*- C++ -*-
#ifndef HERWIG_MEPP2VV_H
#define HERWIG_MEPP2VV_H
//
// This is the declaration of the MEPP2VV class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2VV class implements the production of \f$W^+W^-\f$,
 * \f$W^\pm Z^0\f$ and \f$Z^0Z^o\f$ in hadron-hadron collisions.
 *
 * @see \ref MEPP2VVInterfaces "The interfaces"
 * defined for MEPP2VV.
 */
class MEPP2VV: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2VV();

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
   * Return the process being run (WW/ZZ/WZ).
   */
  virtual int process() const { return process_; }

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
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(double cthmin, double cthmax, const double r);

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
   * Matrix element for \f$f\bar{f}\to W^+W^-\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  The first  outgoing W polarization vectors
   * @param v2  The second outgoing W polarization vectors
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double WWME(vector<SpinorWaveFunction>    & f1,
	      vector<SpinorBarWaveFunction> & a1,
	      vector<VectorWaveFunction>    & v1,
	      vector<VectorWaveFunction>    & v2,
	      bool me) const;
  
  /**
   * Matrix element for \f$f\bar{f}\to W^\pm Z^0\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  The outgoing W polarization vectors
   * @param v2  The outgoing Z polarization vectors
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double WZME(vector<SpinorWaveFunction>    & f1,
	      vector<SpinorBarWaveFunction> & a1,
	      vector<VectorWaveFunction>    & v1,
	      vector<VectorWaveFunction>    & v2,
	      bool me) const;
  
  /**
   * Matrix element for \f$f\bar{f}\to Z^0Z^0\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  The first  outgoing Z polarization vectors 
   * @param v2  The second outgoing Z polarization vectors
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double ZZME(vector<SpinorWaveFunction>    & f1,
	      vector<SpinorBarWaveFunction> & a1,
	      vector<VectorWaveFunction>    & v1,
	      vector<VectorWaveFunction>    & v2,
	      bool me) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VV & operator=(const MEPP2VV &) = delete;

private:

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   *   FFWVertex
   */
  AbstractFFVVertexPtr FFWvertex_;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr FFZvertex_;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr WWWvertex_;
  //@}

  /**
   *  Processes
   */
  unsigned int process_;

  /**
   *  Allowed flavours of the incoming quarks
   */
  int maxflavour_;

  /**
   *  Treatment of the the boson masses
   */
  unsigned int massOption_;

  /**
   *  The matrix element
   */
  mutable ProductionMatrixElement me_;
};

}

#endif /* HERWIG_MEPP2VV_H */
