// -*- C++ -*-
#ifndef Herwig_ReweightEW_H
#define Herwig_ReweightEW_H
//
// This is the declaration of the ReweightEW class.
//

#include "ThePEG/MatrixElement/ReweightBase.h"
#include <boost/array.hpp>

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the ReweightEW class.
 *
 * @see \ref ReweightEWInterfaces "The interfaces"
 * defined for ReweightEW.
 */
class ReweightEW: public ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ReweightEW();

  /**
   * The destructor.
   */
  virtual ~ReweightEW();
  //@}

public:

  /**
   * Return the weight for the kinematical configuation provided by
   * the assigned XComb object (in the LastXCombInfo base class).
   */
  virtual double weight() const;

  /**
   * Return values of the last evaluation (double/doubles in GeV2)
   */
  double lastS() const {return thelasts;} 
  double lastT() const {return thelastt;} 
  double lastK() const {return thelastk;} 

  void setSTK(const double &s, const double &t, const double &K); 


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


private:

  /**
   * The last s
   */
  mutable double thelasts;

  /**
   * The last t
   */
  mutable double thelastt;

  /**
   * The last K-factor
   */
  mutable double thelastk;

  /**
   * The table of K factors to be read from file 
   */

  // tab[40000][5];

  boost::array<boost::array<double,6>,40001> tab;
  //  boost::array<boost::array<double,6>,250001> tab;

  /**
   *  EW K factor filename
   */
  string filename;

public:
  /**
   * Computation of K factors from table (s and t in GeV)
   */
  double EWKFac(unsigned int f, double s, double t) const;

private:
  /**
   * initialize tables
   */  
  void inittable();

private:

  /**
   * Describe a concrete base class with persistent data.
   */
  static ClassDescription<ReweightEW> initReweightEW;

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ReweightEW & operator=(const ReweightEW &);

};

}

#endif /* Herwig_ReweightEW_H */
