// -*- C++ -*-
#ifndef Herwig_MEPPto3P0Jet_H
#define Herwig_MEPPto3P0Jet_H
//
// This is the declaration of the MEPPto3P0Jet class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "OniumParameters.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPPto3P0Jet class implements the colour singlet processes for
 * \f$gq\to^3\!\!P_0 q\f$, \f$q\bar{q}\to^3\!\!P_0 g\f$ and
 * \f$gg\to^3\!\!P_0 g\f$.
 *
 * @see \ref MEPPto3P0JetInterfaces "The interfaces"
 * defined for MEPPto3P0Jet.
 */
class MEPPto3P0Jet: public HwMEBase {

public:
  
  /**
   * The default constructor.
   */
  MEPPto3P0Jet() : O1_(ZERO), state_(ccbar), n_(1), process_(0), mOpt_(0)
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return 3;
  }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const {
    return 0;
  }

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPPto3P0Jet & operator=(const MEPPto3P0Jet &) = delete;

private:
  
  /**
   *  Access to the parameters for the quarkonium states
   */
  OniumParametersPtr params_;
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy5 O1_;

  /**
   *  Type of state
   */
  OniumState state_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;

  /**
   *  Which processes to generate
   */
  unsigned int process_;

  /**
   *  Option for the onium mass
   */
  unsigned int mOpt_;
};

}

#endif /* Herwig_MEPPto3P0Jet_H */
