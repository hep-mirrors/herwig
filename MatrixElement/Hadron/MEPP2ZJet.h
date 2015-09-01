// -*- C++ -*-
#ifndef HERWIG_MEPP2ZJet_H
#define HERWIG_MEPP2ZJet_H
//
// This is the declaration of the MEPP2ZJet class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2ZJet class implements the matrix element for the production
 * of a Z boson + a jet including the decay of the Z including \f$Z/\gamma\f$ 
 * interference
 *
 * @see \ref MEPP2ZJetInterfaces "The interfaces"
 * defined for MEPP2ZJet.
 */
class MEPP2ZJet: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2ZJet();

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

  /**
   *  Matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to Z/\gamma g\f$.
   * @param fin   Spinors for incoming quark
   * @param ain   Spinors for incoming antiquark
   * @param gout  Polarization vectors for the outgoing gluon
   * @param lm    Spinors for outgoing lepton
   * @param lp    Spinors for outgoing antilepton
   * @param me    Whether or not to calculate the matrix element for spin correlations
   **/
  InvEnergy2 qqbarME(vector<SpinorWaveFunction> & fin,
		     vector<SpinorBarWaveFunction> & ain,
		     vector<VectorWaveFunction> & gout,
		     vector<SpinorBarWaveFunction> & lm, vector<SpinorWaveFunction> & lp,
		     bool me=false) const;

  /**
   * Matrix element for \f$qg\to Z/\gamma q\f$.
   * @param fin  Spinors for incoming quark
   * @param gin  Polarization vectors for the incoming gluon
   * @param fout Spinors for outgoing quark
   * @param lm    Spinors for outgoing lepton
   * @param lp    Spinors for outgoing antilepton
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
  InvEnergy2 qgME(vector<SpinorWaveFunction> & fin,vector<VectorWaveFunction> & gin,
		  vector<SpinorBarWaveFunction> & fout,
		  vector<SpinorBarWaveFunction> & lm, vector<SpinorWaveFunction> & lp,
		  bool me=false) const;
  
  /**
   * Matrix element for \f$\bar{q}g\to Z/\gamma\bar{q}\f$.
   * @param fin  Spinors for incoming antiquark
   * @param gin  Polarization vectors for the incoming gluon
   * @param fout Spinors for outgoing antiquark
   * @param lm    Spinors for outgoing lepton
   * @param lp    Spinors for outgoing antilepton
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
  InvEnergy2 qbargME(vector<SpinorBarWaveFunction> & fin,
		     vector<VectorWaveFunction> & gin,
		     vector<SpinorWaveFunction> & fout,
		     vector<SpinorBarWaveFunction> & lm, vector<SpinorWaveFunction> & lp,
		     bool me=false) const;
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2ZJet> initMEPP2ZJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2ZJet & operator=(const MEPP2ZJet &);

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
   *  Pointer to the photon vertex
   */
  AbstractFFVVertexPtr _theFFPVertex;

  /**
   *  Pointer to the \f$qqg\f$ vertex
   */
  AbstractFFVVertexPtr _theQQGVertex;
  //@}

  /**
   *  @name Pointers to the \f$Z^0\f$ and \f$\gamma\f$ ParticleData objects
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr _z0;

  /**
   *  Pointer to the photon ParticleData object
   */
  tcPDPtr _gamma;
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
  int _maxflavour;

  /**
   *  Control over which Z decay modes to include
   */
  int _zdec;

  /**
   *  Which terms to include
   */
  unsigned int _gammaZ;

  /**
   *  Option for the treatment of the W off-shell effects
   */
  unsigned int _widthopt;
  //@}

  /**
   *  Probability of selecting \f$1/s^2\f$ for the jacobian
   *  transformation of the boson mass
   */
  double _pprob;
  
  /**
   * Matrix element for spin correlations
   */
  mutable ProductionMatrixElement _me;

  /**
   *  Storage of the scale to avoid the need to recalculate
   */
  Energy2 _scale;
  
  /**
   *  Storage of the off-shell Z mass to avoid the need to recalculate
   */
  Energy2 _mz2;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2ZJet. */
template <>
struct BaseClassTrait<Herwig::MEPP2ZJet,1> {
  /** Typedef of the first base class of MEPP2ZJet. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2ZJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2ZJet>
  : public ClassTraitsBase<Herwig::MEPP2ZJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2ZJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2ZJet is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2ZJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2ZJet_H */
