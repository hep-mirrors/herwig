// -*- C++ -*-
#ifndef Herwig_RPVSSSVertex_H
#define Herwig_RPVSSSVertex_H
//
// This is the declaration of the RPVSSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"
#include "Herwig/Models/Susy/MixingMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVSSSVertex class.
 *
 * @see \ref RPVSSSVertexInterfaces "The interfaces"
 * defined for RPVSSSVertex.
 */
class RPVSSSVertex: public Helicity::SSSVertex {

public:

  /**
   * The default constructor.
   */
  RPVSSSVertex();

  /**
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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
  RPVSSSVertex & operator=(const RPVSSSVertex &) = delete;

private:

  /**
   *  Which types of interaction to include
   */
  unsigned int interactions_;

  /**
   * The scale at which the coupling was last evaluated.  
   */
  Energy2 q2Last_;
  
  /**
   * The value of \f$ \sqrt{4\pi\alpha}\f$ when it was last evaluated.
   */
  double gLast_;

  /**
   *  Sfermion mixing
   */
  //@{
  /**
   *  Stau mixing
   */
  tMixingMatrixPtr stau_;

  /**
   *  Sbottom mixing
   */
  tMixingMatrixPtr sbottom_;

  /**
   *  Stop mixing
   */
  tMixingMatrixPtr stop_;
  //@}

  /**
   *  Couplings of the scalar Higgs bosons to other scalars
   */
  vector<vector<vector<complex<Energy> > > > scalarScalarScalar_;

  /**
   *  Couplings of the scalar Higgs bosons to pseudoscalars
   */
  vector<vector<vector<complex<Energy> > > > scalarPseudoPseudo_;

  /**
   *  Couplings of the neutral scalar Higgs bosons to two charged Higgs bosons
   */
  vector<vector<vector<complex<Energy> > > > scalarChargedCharged_;

  /**
   *  Couplings of the neutral pseudoscalar Higgs bosons to two charged Higgs bosons
   */
  vector<vector<vector<complex<Energy> > > > pseudoChargedCharged_;

  /**
   *  Couplings of the scalar Higgs bosons to up squarks
   */
  vector<vector<vector<vector<complex<Energy> > > > > scalarSup_;

  /**
   *  Couplings of the scalar Higgs bosons to down squarks
   */
  vector<vector<vector<vector<complex<Energy> > > > > scalarSdown_;

  /**
   *  Couplings of the pseudo Higgs bosons to up squarks
   */
  vector<vector<complex<Energy> > > pseudoSup_;

  /**
   *  Couplings of the pseudo Higgs bosons to down squarks
   */
  vector<vector<complex<Energy> > > pseudoSdown_;

  /**
   *  Couplings of the scalar Higgs bosons to sneutinos
   */
  vector<vector<complex<Energy> > > scalarSneutrino_;

  /**
   *  Couplings of the scalar Higgs bosons to charged sleptons
   */
  vector<vector<vector<vector<complex<Energy> > > > > scalarSlepton_;

  /**
   *  Couplings of the pseudo Higgs bosons to charged sleptons
   */
  vector<vector<complex<Energy> > > pseudoSlepton_;

  /**
   *  Couplings of the charged Higgs to squarks
   */
  vector<vector<vector<vector<complex<Energy> > > > > chargedSquark_;

  /**
   *  Couplings of the charged Higgs to sleptons
   */
  vector<vector<vector<complex<Energy> > > > chargedSlepton_;

};

}

#endif /* Herwig_RPVSSSVertex_H */
