// -*- C++ -*-
#ifndef HERWIG_MEqq2gZ2ffNason_H
#define HERWIG_MEqq2gZ2ffNason_H
//
// This is the declaration of the MEqq2gZ2ffNason class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.fh"
#include "Herwig++/Utilities/Statistic.h"
#include "MEqq2gZ2ffNason.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEqq2gZ2ffNason class implements the products of Standard Model
 * fermion antifermion pairs via the \f$Z^0\f$ resonance including
 * photon interference terms.
 *
 * @see \ref MEqq2gZ2ffNasonInterfaces "The interfaces"
 * defined for MEqq2gZ2ffNason.
 */
class MEqq2gZ2ffNason: public ME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEqq2gZ2ffNason();

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
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

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

  /**
   * Matrix element for \f$q\bar{q}\to \gamma/Z \to f\bar{f}\f$.
   * @param fin  Spinors for incoming quark
   * @param ain  Spinors for incoming antiquark
   * @param fout Spinors for incoming quark
   * @param aout Spinors for incoming antiquark
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction>    & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<SpinorBarWaveFunction> & fout,
		 vector<SpinorWaveFunction>    & aout,
		 bool me) const;

  /**
   *  Calculate the correction weight
   */
  double NLOweight() const;

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
  virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEqq2gZ2ffNason> initMEqq2gZ2ffNason;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2gZ2ffNason & operator=(const MEqq2gZ2ffNason &);

private:

  /**
   *  Pointer to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the Z vertex
   */
  FFVVertexPtr _theFFZVertex;

  /**
   *  Pointer to the photon vertex
   */
  FFVVertexPtr _theFFPVertex;
  //@}

  /**
   *  Pointers to the intermediate resonances
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
   *  Switches to control the particles in the hard process
   */
  //@{
  /**
   *  Allowed flavours for the incoming quarks
   */
  unsigned int _maxflavour;

  /**
   *  Whether to include both \f$Z^0\f$ and \f$\gamma\f$ or only one
   */
  unsigned int _gammaZ;
  
  /**
   *  Which processes to include
   */
  unsigned int _process;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$x\f$ variable
   */
  double _x;

  /**
   *  The \f$\tilde{v}
   */
  double _v;
  //@}

  /**
   * Statistics on the weights for testing
   */
  mutable vector<Statistic> _posx,_negx,_posv,_negv;
  mutable vector<Statistic> _posxp,_negxp,_posvp,_negvp;
  mutable vector<Statistic> _posxn,_negxn,_posvn,_negvn;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2gZ2ffNason. */
template <>
struct BaseClassTrait<Herwig::MEqq2gZ2ffNason,1> {
  /** Typedef of the first base class of MEqq2gZ2ffNason. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2gZ2ffNason class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEqq2gZ2ffNason>
  : public ClassTraitsBase<Herwig::MEqq2gZ2ffNason> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEqq2gZ2ffNason"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MEqq2gZ2ffNason class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwNasonME.so"; }
};

/** @endcond */

}

#include "MEqq2gZ2ffNason.icc"

#endif /* HERWIG_MEqq2gZ2ffNason_H */
