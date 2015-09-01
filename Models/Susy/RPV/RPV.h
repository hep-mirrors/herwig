// -*- C++ -*-
#ifndef HERWIG_RPV_H
#define HERWIG_RPV_H
//
// This is the declaration of the RPV class.
//

#include "Herwig/Models/Susy/MSSM.h"
#include "RPV.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPV class.
 *
 * @see \ref RPVInterfaces "The interfaces"
 * defined for RPV.
 */
class RPV: public MSSM {

public:

  /**
   * The default constructor.
   */
  RPV() : lambdaLLE_(3,vector<vector<double> >(3,vector<double>(3,0.))),
	  lambdaLQD_(3,vector<vector<double> >(3,vector<double>(3,0.))),
	  lambdaUDD_(3,vector<vector<double> >(3,vector<double>(3,0.))),
	  triLinearOnly_(false)
  {}


  /**
   *  LLE couplings
   */
  const vector<vector<vector<double> > > & lambdaLLE() {return lambdaLLE_;}

  /**
   * LQD couplings
   */
  const vector<vector<vector<double> > > & lambdaLQD() {return lambdaLQD_;}

  /**
   * UDD couplings
   */
  const vector<vector<vector<double> > > & lambdaUDD() {return lambdaUDD_;}

  /**
   *  Sneutrino vevs
   */
  const vector<Energy> & sneutrinoVEVs() {return vnu_;}

  /**
   *  Bilinear coupings
   */
  const vector<Energy> & epsilon() {return epsilon_;}

  /**
   *  Blinear soft terms
   */
  const vector<Energy> & epsilonB() {return epsB_;}

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
  virtual void doinit();
  //@}

protected:

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RPV & operator=(const RPV &);

private:

  /**
   *  Trillinear couplings
   */
  //@
  /**
   *  LLE couplings
   */
  vector<vector<vector<double> > > lambdaLLE_;

  /**
   * LQD couplings
   */
  vector<vector<vector<double> > > lambdaLQD_;

  /**
   * UDD couplings
   */
  vector<vector<vector<double> > > lambdaUDD_;
  //@

  /**
   *  Sneutrino vevs
   */
  vector<Energy> vnu_;

  /**
   *  The bilinear parameters
   */
  vector<Energy> epsilon_;

  /**
   *  The bilinear soft terms
   */
  vector<Energy> epsB_;

  /**
   *  Squark mixing matrices
   */
  //@{
  /**
   *  For up-type squarks
   */
  MixingMatrixPtr upSquarkMix_;

  /**
   *  For down-type squarks
   */
  MixingMatrixPtr downSquarkMix_;
  //@}

  /**
   *  New vertices
   */
  //@{
  /**
   *  LLE vertex
   */
  AbstractFFSVertexPtr LLEVertex_;

  /**
   *  LQD vertex
   */
  AbstractFFSVertexPtr LQDVertex_;

  /**
   *  UDD vertex
   */
  AbstractFFSVertexPtr UDDVertex_;
  //@}

  /**
   *  Switch to only do trilinears
   */
  bool triLinearOnly_;
};

}

#endif /* HERWIG_RPV_H */
