// -*- C++ -*-
//
// MEPP2GammaGamma.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2GammaGamma_H
#define HERWIG_MEPP2GammaGamma_H
//
// This is the declaration of the MEPP2GammaGamma class.
//

#include "Herwig++/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2GammaGamma class implements the production of photon pairs in 
 * hadron hadron collisions.
 *
 * @see \ref MEPP2GammaGammaInterfaces "The interfaces"
 * defined for MEPP2GammaGamma.
 */
class MEPP2GammaGamma: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2GammaGamma() : _maxflavour(5),_process(0) {
    massOption(vector<unsigned int>(2,0));
  }
  
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
   *  Members to return the matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to \gamma\gamma\f$.
   * @param fin Spinors for incoming quark
   * @param ain Spinors for incoming antiquark
   * @param p1  Polarization vectors for the first  outgoing photon
   * @param p2  Polarization vectors for the second outgoing photon
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction> & fin, vector<SpinorBarWaveFunction> & ain,
		 vector<VectorWaveFunction> & p1 , vector<VectorWaveFunction>    & p2 ,
		 bool me) const;

  /**
   * Matrix element for \f$gg \to \gamma\gamma\f$.
   * @param g1  Polarization vectors for the first  incoming gluon
   * @param g2  Polarization vectors for the second incoming gluon
   * @param p1  Polarization vectors for the first  outgoing photon
   * @param p2  Polarization vectors for the second outgoing photon
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double ggME(vector<VectorWaveFunction> & g1 , vector<VectorWaveFunction>    & g2 ,
	      vector<VectorWaveFunction> & p1 , vector<VectorWaveFunction>    & p2 ,
	      bool me) const;
  //@}

  /**
   *  \f$gg\to\gamma\gamma\f$ matrix element for the \f$++++\f$ helicity configuration.
   * @param s The \f$s\f$ invariant
   * @param t The \f$t\f$ invariant
   * @param u The \f$u\f$ invariant
   */
  Complex ggme(Energy2 s,Energy2 t,Energy2 u) const {
    double ltu(log(abs(t/u)));
    double frac1((t-u)/s),frac2((sqr(t)+sqr(u))/sqr(s));
    double thetatu = (t/u<0) ? 0 : 1;
    double thetat  = (t<ZERO)   ? 0 : 1;
    double thetau  = (u<ZERO)   ? 0 : 1;
    using Constants::pi;
    return Complex(1.+frac1*ltu+0.5*frac2*(sqr(ltu)+sqr(pi)*thetatu),
		   -pi*(thetat-thetau)*(frac1+frac2*ltu));
  }

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2GammaGamma> initMEPP2GammaGamma;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2GammaGamma & operator=(const MEPP2GammaGamma &);

private:

  /**
   *  Pointer to the quark-antiquark-photon vertex
   */
  AbstractFFVVertexPtr _photonvertex;

  /**
   *  Maximum PDG code of the quarks allowed
   */
  unsigned int _maxflavour;

  /**
   *  Option for which processes to include
   */
  unsigned int _process;
  
  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   *  weights for the different quark annhilation diagrams
   */
  mutable double _diagwgt[2];

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2GammaGamma. */
template <>
struct BaseClassTrait<Herwig::MEPP2GammaGamma,1> {
  /** Typedef of the first base class of MEPP2GammaGamma. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2GammaGamma class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2GammaGamma>
  : public ClassTraitsBase<Herwig::MEPP2GammaGamma> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2GammaGamma"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2GammaGamma is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2GammaGamma depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2GammaGamma_H */
