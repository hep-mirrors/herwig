// -*- C++ -*-
#ifndef HERWIG_SextetModel_H
#define HERWIG_SextetModel_H
//
// This is the declaration of the SextetModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
#include "SextetModel.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Models
 * 
 *  This class is used instead of the StandardModel class for the 
 *
 *
 * @see \ref SextetModelInterfaces "The interfaces"
 * defined for SextetModel.
 */
class SextetModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  SextetModel() : g1L_(3,0.), g1R_(3,0.), g1pR_(3,0.), g1ppR_(3,0.),
		  g2_(3,0.), g2p_(3,0.), g3L_(3,0.),
		  enableScalarSingletY43_(false),enableScalarSingletY13_(false),
		  enableScalarSingletY23_(false),enableScalarTripletY13_(false),
		  enableVectorDoubletY16_(false),enableVectorDoubletY56_(false) {
    useMe();
  }

  /**
   *  Access to the couplings
   */
  //@{
  /**
   * The \f$SU(2)\f$ quark-doublet coupling to \f$\Phi_{6,1,1/3}\f$
   */
  const vector<double> & g1L() const {return g1L_;}

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,1/3}\f$
   */
  const vector<double> & g1R() const {return g1R_;}

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,-2/3}\f$
   */
  const vector<double> & g1pR() const {return g1pR_;}

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,4/3}\f$
   */
  const vector<double> & g1ppR() const {return g1ppR_;}

  /**
   * The coupling to \f$V^\mu_{6,2,-1/6}\f$
   */
  const vector<double> & g2() const {return g2_;}

  /**
   * The coupling to \f$V^\mu_{6,2,5/6}\f$
   */
  const vector<double> & g2p() const {return g2p_;}

  /**
   * Coupling to \f$\Phi_{6,3,1/3}\f$
   */
  const vector<double> & g3L() const {return g3L_;}
  //@}

  /**
   *  Switches to decide which particles to include
   */
  //@{
  /**
   * Scalar Singlet \f$Y = 4/3\f$
   */
  bool ScalarSingletY43Enabled() const {return enableScalarSingletY43_;}

  /**
   * Scalar Singlet \f$Y = -1/3\f$
   */
  bool ScalarSingletY13Enabled() const {return enableScalarSingletY13_;}

  /**
   * Scalar Singlet \f$Y = -2/3\f$
   */
  bool ScalarSingletY23Enabled() const {return enableScalarSingletY23_;}

  /**
   * Scalar Triplet \f$Y =  1/3\f$
   */
  bool ScalarTripletY13Enabled() const {return enableScalarTripletY13_;}

  /**
   * Vector Doublet \f$Y = -1/6\f$
   */
  bool VectorDoubletY16Enabled() const {return enableVectorDoubletY16_;}

  /**
   * Vector Doublet \f$Y =  5/6\f$
   */
  bool VectorDoubletY56Enabled() const {return enableVectorDoubletY56_;}
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

  /**
   *  Member to implement the command to enable particular diquarks
   */
  string doEnable(string command);
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SextetModel & operator=(const SextetModel &);

private:

  /**
   *  Pointers to the vertex objects
   */
  //@{
  /**
   *  Pointer to the object handling the strong coupling of a 
   *  vector sextet to one gluon
   */
  AbstractVVVVertexPtr VVVVertex_;

  /**
   *  Pointer to the object handling the strong coupling of a 
   *  vector sextet to two gluons
   */
  AbstractVVVVVertexPtr VVVVVertex_;

  /**
   *  Pointer to the object handling the strong coupling of a 
   *  scalar sextet to one gluon
   */
  AbstractVSSVertexPtr VSSVertex_;

  /**
   *  Pointer to the object handling the strong coupling of a 
   *  scalar sextet to two gluons
   */
  AbstractVVSSVertexPtr VVSSVertex_;

  /**
   *  Pointer to the object handling the coupling of two quarks
   *  to a vector sextet
   */
  AbstractFFVVertexPtr FFVVertex_;

  /**
   *  Pointer to the object handling the coupling of two quarks
   *  to a scalar sextet
   */
  AbstractFFSVertexPtr FFSVertex_;
  //@}

  /**
   *  Couplings
   */
  //@{
  /**
   * The \f$SU(2)\f$ quark-doublet coupling to \f$\Phi_{6,1,1/3}\f$
   */
  vector<double> g1L_;

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,1/3}\f$
   */
  vector<double> g1R_;

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,-2/3}\f$
   */
  vector<double> g1pR_;

  /**
   * The \f$SU(2)\f$ singlet coupling to \f$\Phi_{6,1,4/3}\f$
   */
  vector<double> g1ppR_;

  /**
   * The coupling to \f$V^\mu_{6,2,-1/6}\f$
   */
  vector<double> g2_;

  /**
   * The coupling to \f$V^\mu_{6,2,5/6}\f$
   */
  vector<double> g2p_;

  /**
   * Coupling to \f$\Phi_{6,3,1/3}\f$
   */
  vector<double> g3L_;
  //@}

  /**
   *  Switches to decide which particles to include
   */
  //@{
  /**
   * Scalar Singlet \f$Y = 4/3\f$
   */
  bool enableScalarSingletY43_;

  /**
   * Scalar Singlet \f$Y = -1/3\f$
   */
  bool enableScalarSingletY13_;

  /**
   * Scalar Singlet \f$Y = -2/3\f$
   */
  bool enableScalarSingletY23_;

  /**
   * Scalar Triplet \f$Y =  1/3\f$
   */
  bool enableScalarTripletY13_;

  /**
   * Vector Doublet \f$Y = -1/6\f$
   */
  bool enableVectorDoubletY16_;

  /**
   * Vector Doublet \f$Y =  5/6\f$
   */
  bool enableVectorDoubletY56_;
  //@}
};

}

#endif /* HERWIG_SextetModel_H */
