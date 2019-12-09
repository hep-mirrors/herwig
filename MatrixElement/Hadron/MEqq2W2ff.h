// -*- C++ -*-
//
// MEqq2W2ff.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEqq2W2ff_H
#define HERWIG_MEqq2W2ff_H
//
// This is the declaration of the MEqq2W2ff class.
//

#include "Herwig/MatrixElement/DrellYanBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEqq2W2ff class implements the matrix element for \f$q\bar{q'}\to W^\pm\f$
 * including the decay of the \f$W^\pm\f$ to Standard Model fermions.
 *
 * @see \ref MEqq2W2ffInterfaces "The interfaces"
 * defined for MEqq2W2ff.
 */
class MEqq2W2ff: public DrellYanBase {

public:

  /**
   * The default constructor.
   */
  MEqq2W2ff();

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2W2ff & operator=(const MEqq2W2ff &) = delete;

private:

  /**
   *  Pointer to the W vertex
   */
  AbstractFFVVertexPtr  _theFFWVertex;

  /**
   *  Pointers to the intermediates resonances
   */
  //@{
  /**
   *  Pointer to the \f$W^+\f$
   */
  tcPDPtr _wp;

  /**
   *  Pointer to the \f$W^-\f$
   */
  tcPDPtr _wm;
  //@}

  /**
   *  Switches to control the particles in the hard process
   */
  //@{
  /**
   *  The allowed flavours of the incoming quarks
   */
  unsigned int _maxflavour;

  /**
   *  Which intermediate \f$W^\pm\f$ bosons to include
   */
  unsigned int _plusminus;

  /**
   *  Which decay products of the \f$W^\pm\f$ to include
   */
  unsigned int _process;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;
};

}

#endif /* HERWIG_MEqq2W2ff_H */
