// -*- C++ -*-
//
// MEPP2GammaJet.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2GammaJet_H
#define HERWIG_MEPP2GammaJet_H
//
// This is the declaration of the MEPP2GammaJet class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup MatrixElements
 * The MEPP2GammaJet class implements the matrix element for photon+jet 
 * production in hadron-hadron collisions.
 *
 * @see \ref MEPP2GammaJetInterfaces "The interfaces"
 * defined for MEPP2GammaJet.
 */
class MEPP2GammaJet: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2GammaJet();

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

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

private:

  /**
   *  Members to return the matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to g\gamma\f$.
   * @param fin   Spinors for incoming quark
   * @param ain   Spinors for incoming antiquark
   * @param gout  Polarization vectors for the outgoing gluon
   * @param pout  Polarization vectors for the outgoing photon
   * @param me    Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction> & fin, vector<SpinorBarWaveFunction> & ain,
		 vector<VectorWaveFunction> & gout, vector<VectorWaveFunction>    & pout,
		 bool me) const;

  /**
   * Matrix element for \f$qg\to \gamma q\f$.
   * @param fin   Spinors for incoming quark
   * @param gin   Polarization vectors for the incoming gluon
   * @param pout  Polarization vectors for the outgoing photon
   * @param fout  Spinors for outgoing quark
   * @param me    Whether or not to calculate the matrix element for spin correlations
   */
  double qgME(vector<SpinorWaveFunction> & fin,vector<VectorWaveFunction>     & gin,
	      vector<VectorWaveFunction> & pout,vector<SpinorBarWaveFunction> & fout,
	      bool me) const;

  /**
   * Matrix element for \f$\bar{q}g\to \gamma \bar{q}\f$.
   * @param ain   Spinors for the incoming antiquark
   * @param gin   Polarization vectors for the incoming gluon
   * @param pout  Polarization vectors for the outgoing photon
   * @param aout  Spinors for the outgoing antiquark
   * @param me    Whether or not to calculate the matrix element for spin correlations
   */
  double qbargME(vector<SpinorBarWaveFunction> & ain, vector<VectorWaveFunction> & gin,
		 vector<VectorWaveFunction> & pout, vector<SpinorWaveFunction> & aout,
		 bool me) const;
  //@}
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2GammaJet & operator=(const MEPP2GammaJet &) = delete;

private:

  /**
   *  Pointer to the quark-antiquark-gluon vertex
   */
  AbstractFFVVertexPtr _gluonvertex;

  /**
   *  Pointer to the quark-antiquark-photon vertex
   */
  AbstractFFVVertexPtr _photonvertex;

  /**
   *  Maximum PDG code of the quarks allowed
   */
  int _maxflavour;

  /**
   *  Option for which processes to include
   */
  unsigned int _processopt;
  
  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

  /**
   *  Scale prefactor
   */
  double scalePreFactor_;
};

}

#endif /* HERWIG_MEPP2GammaJet_H */
