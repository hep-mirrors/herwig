// -*- C++ -*-
#ifndef Herwig_MESimpleZJet_H
#define Herwig_MESimpleZJet_H
//
// This is the declaration of the MESimpleZJet class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MESimpleZJet class.
 *
 * @see \ref MESimpleZJetInterfaces "The interfaces"
 * defined for MESimpleZJet.
 */
class MESimpleZJet: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MESimpleZJet();

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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

  /**
   *  Matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to Z^0 g\f$.
   * @param fin   Spinors for incoming quark
   * @param ain   Spinors for incoming antiquark
   * @param gout  Polarization vectors for the outgoing gluon
   * @param Zout  Polarization vectors for the outgoing Z
   * @param me    Whether or not to calculate the matrix element for spin correlations
   **/
  double qqbarME(std::vector<SpinorWaveFunction> & fin,
		 std::vector<SpinorBarWaveFunction> & ain,
		 std::vector<VectorWaveFunction> & gout,
		 std::vector<VectorWaveFunction> & Zout,
		 bool me=false) const;

  /**
   * Matrix element for \f$qg\to Z^0 q\f$.
   * @param fin  Spinors for incoming quark
   * @param gin  Polarization vectors for the incoming gluon
   * @param fout Spinors for outgoing quark
   * @param Zout Polarization vectors for the outgoing Z
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
   double qgME(std::vector<SpinorWaveFunction> & fin,
             std::vector<VectorWaveFunction> & gin,
             std::vector<SpinorBarWaveFunction> & fout,
             std::vector<VectorWaveFunction> & Zout,
             bool me=false) const;

  /**
   * Matrix element for \f$\bar{q}g\to Z^0\bar{q}\f$.
   * @param fin  Spinors for incoming antiquark
   * @param gin  Polarization vectors for the incoming gluon
   * @param fout Spinors for outgoing antiquark
   * @param Zout Polarization vectors for the outgoing Z
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
   double qbargME(std::vector<SpinorBarWaveFunction> & fin,
                std::vector<VectorWaveFunction> & gin,
                std::vector<SpinorWaveFunction> & fout,
                std::vector<VectorWaveFunction> & Zout,
                bool me=false) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MESimpleZJet & operator=(const MESimpleZJet &);

private:

  /**
   *  Vertices for the helicity amplitude calculation
   */
  //@{
  /**
   *  Pointer to the Z vertex
   */
  AbstractFFVVertexPtr _theFFZVertex;

  /**
   *  Pointer to the \f$qqg\f$ vertex
   */
  AbstractFFVVertexPtr _theQQGVertex;
  //@}

  /**
   *  @name Pointers to the Z ParticleData objects
   */
  //@{
  /**
   *  The \f$Z^0\f$ data pointer
   */
  tcPDPtr _z0;
  //@}

  /**
   * @name Switches to control the particles in the hard process
   */
  //@{
  /**
   *  Subprocesses to include
   */
  unsigned int _process;

  /**
   *  Allowed flavours for the incoming quarks
   */
  unsigned int _maxflavour;

  /**
   * Matrix element for spin correlations
   */
  mutable ProductionMatrixElement _me;

};

}

#endif /* Herwig_MESimpleZJet_H */
