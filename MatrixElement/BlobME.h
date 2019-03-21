// -*- C++ -*-
#ifndef Herwig_BlobME_H
#define Herwig_BlobME_H
//
// This is the declaration of the BlobME class.
//

#include "ThePEG/MatrixElement/BlobMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \brief BlobME serves as a base class for special processes such as
 * instanton or sphaleron induced ones.
 *
 * \author Simon Platzer
 *
 * @see \ref BlobMEInterfaces "The interfaces"
 * defined for BlobME.
 */
class BlobME: public BlobMEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BlobME();

  /**
   * The destructor.
   */
  virtual ~BlobME();
  //@}

public:

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const {
    return lastSHat();
  }

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector.
   */
  virtual bool generateKinematics(const double * r) {
    jacobian(thePhasespace->generateTwoToNKinematics(r,meMomenta()));
    return jacobian() > 0.0;
  }
  
  /**
   * Return the additional number of objects to generate
   */
  size_t nAdditional() const { return theNAdditional; }

  /**
   * Set the additional number of objects to generate
   */
  void nAdditional(size_t n) { theNAdditional = n; }

  /**
   * Return the minimal number of final state particles
   */
  virtual size_t nOutgoing() const = 0;

  /**
   * The number of internal degreed of freedom used in the matrix
   * element. This default version returns 0;
   */
  virtual int nDim() const {
    return thePhasespace->nDimPhasespace(nOutgoing()+nAdditional());
  }

  /**
   * Set the xcomb object
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Return the phase space used
   */
  Ptr<MatchboxPhasespace>::tptr phasespace() const { return thePhasespace; }

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The phase space to be used
   */
  Ptr<MatchboxPhasespace>::ptr thePhasespace;

  /**
   * The multiplicity of additional objects to consider
   */
  size_t theNAdditional;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BlobME & operator=(const BlobME &) = delete;

};

}

#endif /* Herwig_BlobME_H */
