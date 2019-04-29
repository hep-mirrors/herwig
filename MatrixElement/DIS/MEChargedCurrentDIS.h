// -*- C++ -*-
#ifndef HERWIG_MEChargedCurrentDIS_H
#define HERWIG_MEChargedCurrentDIS_H
//
// This is the declaration of the MEChargedCurrentDIS class.
//

#include "DISBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEChargedCurrentDIS class provides the matrix elements for
 * charged current DIS.
 *
 *  By default both the incoming and outgong quarks are assumed to be massless
 *  although the mass of the outgoing quark can be included if required. This
 *  option should be used if top production is included.
 *
 * @see \ref MEChargedCurrentDISInterfaces "The interfaces"
 * defined for MEChargedCurrentDIS.
 */
class MEChargedCurrentDIS: public DISBase {

public:

  /**
   * The default constructor.
   */
  MEChargedCurrentDIS();

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
   * Matrix element for \f$\ell q\to \gamma/Z \to \ell q\f$.
   * @param f1 Fermion on lepton line
   * @param a1 Anti-fermion on lepton line
   * @param f2 Fermion on quark line
   * @param a2 Anti-fermion on quark line
   * @param lorder The order of particles on the lepton line
   * @param qorder The order of particles on the quark line
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double helicityME(vector<SpinorWaveFunction>    & f1 ,
		    vector<SpinorWaveFunction>    & f2,
		    vector<SpinorBarWaveFunction> & a1 ,
		    vector<SpinorBarWaveFunction> & a2,
		    bool lorder, bool qorder,
		    bool me) const;

  /**
   *  Calculate the coefficient A for the correlations in the hard
   *  radiation
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  MEChargedCurrentDIS & operator=(const MEChargedCurrentDIS &) = delete;

private:

  /**
   *  Pointer to the vertex for the helicity calculations
   */
  AbstractFFVVertexPtr _theFFWVertex;

  /**
   *  The allowed flavours of the incoming quarks
   */
  unsigned int _maxflavour;

  /**
   *  Option for the mass of the outgoing quarks
   */
  unsigned int _massopt;

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

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
};

}

#endif /* HERWIG_MEChargedCurrentDIS_H */
