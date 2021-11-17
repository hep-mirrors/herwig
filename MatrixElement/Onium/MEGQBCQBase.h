// -*- C++ -*-
#ifndef Herwig_MEGQBCQBase_H
#define Herwig_MEGQBCQBase_H
//
// This is the declaration of the MEGQBCQBase class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "OniumParameters.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEGQBCQBase class provides a base class for the simulation of \f$gc\to B_c b\f$ processes, handling the kinematics.
 *
 * @see \ref MEGQBCQBaseInterfaces "The interfaces"
 * defined for MEGQBCQBase.
 */
class MEGQBCQBase: public HwMEBase {

public:
  
  /**
   * The default constructor.
   */
  MEGQBCQBase(long pid=531) : id_(pid), n_(1), mOpt_(1)
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
  virtual unsigned int orderInAlphaEW() const  {
    return 0;
  }

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

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
  //@}
  
protected:

  /**
   *  Access to principal quantum number
   */
  unsigned int principleQuantumNumber() const {
    return n_;
  }

  /**
   *  Access to the parameters
   */
  OniumParametersPtr oniumParameters() {
    return params_;
  }

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

public:

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEGQBCQBase & operator=(const MEGQBCQBase &) = delete;

private:

  /**
   *    PDG code for the \f$B_c\f$ state
   */
  long id_;

  /**
   *   Principle quantum number for the state
   */
  unsigned int n_;

  /**
   *  Option for the mass handling
   */
  unsigned int mOpt_;

  /**
   *  The quarkonium parameters
   */
  OniumParametersPtr params_;

  /**
   *  Particle data object for the \f$B_c\f$ state
   */
  PDPtr state_;

  /**
   *  Mass generator for the state
   */
  GenericMassGeneratorPtr massGen_;
};

}

#endif /* Herwig_MEGQBCQBase_H */
