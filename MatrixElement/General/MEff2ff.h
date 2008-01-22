// -*- C++ -*-
//
// MEff2ff.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2ff_H
#define HERWIG_MEff2ff_H
//
// This is the declaration of the MEff2ff class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "MEff2ff.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;

/**
 * This is the implementation of the \f$ 2\to 2\f$ matrix element for
 * a \f$ \Psi \Psi \to \Psi \Psi\f$ process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions.
 *
 * @see \ref MEff2ffInterfaces "The Interfaces"
 * defined for MEff2ff.
 * @see GeneralHardME
 */
class MEff2ff: public GeneralHardME {

public:
  
  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

public:

  /**
   * The default constructor.
   */
  inline MEff2ff();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
  
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

private:
  
  /** @name Functions to compute the ProductionMatrixElement. */
  //@{
  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\Psi\bar{\Psi}\f$
   * @param fin Spinors for first incoming particle
   * @param fbin SpinorBar Wavefunctions for second incoming particle
   * @param fbout SpinorBar Wavefunctions for outgoing particle
   * @param fout Spinors for first outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2ffbHeME(SpinorVector & fin, SpinorBarVector & fbin,
	      SpinorBarVector & fbout, SpinorVector & fout,
	      double & me2) const;

  /**
   * Compute the matrix element for \f$\Psi\Psi\to\Psi\Psi\f$
   * @param fin Spinors for first incoming particle
   * @param fin2 Spinors  for second incoming particle
   * @param fbout SpinorBar for first outgoing particle
   * @param fbout2 SpinorBar Wavefunctions for outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ff2ffHeME(SpinorVector & fin, SpinorVector & fin2,
	    SpinorBarVector & fbout, SpinorBarVector & fbout2,
	    double & me2) const;
  
  /**
   * Compute the matrix element for 
   * \f$\bar{\Psi}\bar{\Psi}\to\bar{\Psi}\bar{\Psi}\f$
   * @param fbin SpinorBars for first incoming particle
   * @param fbin2 SpinorBars  for second incoming particle
   * @param fout Spinors for first outgoing particle
   * @param fout2 Spinors Wavefunctions for outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  fbfb2fbfbHeME(SpinorBarVector & fbin, SpinorBarVector & fbin2,
		SpinorVector & fout, SpinorVector & fout2,
		double & me2) const;

  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\lambda\lambda\f$
   * @param fin Spinors for first incoming particle
   * @param fbin SpinorBar Wavefunctions for second incoming particle
   * @param fbout SpinorBar Wavefunctions for first outgoing particle
   * @param fout Spinors for second outgoing particle
   * @param fout2 Spinor Wavefunctions for first outgoing particle
   * @param fbout2 SpinorBar Wavefunctions for second outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement 
  ffb2mfmfHeME(SpinorVector & fin, SpinorBarVector & fbin, 
	       SpinorBarVector & fbout, SpinorVector & fout,
	       SpinorVector & fout2, SpinorBarVector & fbout2,
	       double & me2) const;
  //@}

  /**
   * Construct the vertex information for the spin correlations
   * @param sub Pointer to the relevent SubProcess
   */
  virtual void constructVertex(tSubProPtr sub);

protected:
  
  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;
  
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
  virtual void doinit() throw(InitException);
  //@}

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEff2ff> initMEff2ff;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2ff & operator=(const MEff2ff &);

private:
  
  /**
   * Store the vector of FFSVertex pairs
   */
  vector<pair<AbstractFFSVertexPtr, AbstractFFSVertexPtr> > theScaV;

  /**
   * Store the vector of FFVVertex pairs
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > theVecV;

  /**
   * Store the vector of FFTVertex pairs
   */
  vector<pair<AbstractFFTVertexPtr, AbstractFFTVertexPtr> > theTenV;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2ff. */
template <>
struct BaseClassTrait<Herwig::MEff2ff,1> {
  /** Typedef of the first base class of MEff2ff. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2ff class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2ff>
  : public ClassTraitsBase<Herwig::MEff2ff> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEff2ff"; }
};

/** @endcond */

}

#include "MEff2ff.icc"

#endif /* HERWIG_MEff2ff_H */
