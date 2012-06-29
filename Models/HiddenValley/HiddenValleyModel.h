// -*- C++ -*-
#ifndef HERWIG_HiddenValleyModel_H
#define HERWIG_HiddenValleyModel_H
//
// This is the declaration of the HiddenValleyModel class.
//

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenValleyModel class.
 *
 * @see \ref HiddenValleyModelInterfaces "The interfaces"
 * defined for HiddenValleyModel.
 */
class HiddenValleyModel: public StandardModel {

/**
 *  Enumeration to define the type of the confining group 
 */
  enum GroupType { SU, SO, SP };

public:

  /**
   * The default constructor.
   */
  HiddenValleyModel();

public:

  /**
   *  Properties of the new confining gauge group
   */
  //@{
  /**
   *  The hidden colour charge of the fundamental representation
   */
  double CF() const {
    switch (groupType_) {
    case SU:
      return 0.5*(sqr(double(Nc_))-1.)/double(Nc_);
    case SO:
      if(Nc_==3) {
	return     (Nc_-1);
      }
      else {
	return 0.5*(Nc_-1);
      }
    case SP:
      return 0.25*(Nc_+1);
    default:
      assert(false);
    }
    return 0.;
  }

  /**
   *  The hidden colour charge of the adjoint representation
   */
  double CA() const {
    switch (groupType_) {
    case SU:
      return double(Nc_);
    case SO:
      if(Nc_==3) {
	return  2.*(Nc_-2.);
      }
      else {
	return     (Nc_-2.);
      }
    case SP:
      return   0.5*(Nc_+2.);
    default:
      assert(false);
    }
    return 0.;
  }

  /**
   *  The \f$T_R\f$ colour factor
   */
  double TR() const {
    switch (groupType_) {
    case SU: case SP:
      return 0.5;
      break;
    case SO:
      if(Nc_==3) {
	return  2.;
      }
      else {
	return  1.;
      }
    default:
      assert(false);
    }
    return 0.;
  }
  
  /**
   *  The type of the group
   */ 
  GroupType groupType() const {return groupType_;}

  /**
   * The number of colours
   */
  unsigned int NC() const {return Nc_;}

  /**
   *  The number of fermions charged under the new group
   */
  unsigned int NF() const {return Nf_;}
  //@}

  /**
   *  Properties of the new U(1) group
   */
  //@{
  /**
   *  Charge of the left-handed quarks
   */
  double qL() const {return qL_;}
  
  /**
   *  Charge of the right-handed up-type quarks
   */
  double uR() const {return uR_;}
  
  /**
   *  Charge of the right-handed down-type quarks
   */
  double dR() const {return dR_;}
  
  /**
   *  Charge of the left-handed leptons
   */
  double lL() const {return lL_;}

  /**
   *  Charge of the right-handed charged leptons
   */
  double lR() const {return lR_;}

  /**
   *  Coupling
   */
  double gPrime() const {return gPrime_;}

  /**
   *  The couplings of the quirks
   */
  vector<double> qCharge() const {return qCharge_;}
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HiddenValleyModel> initHiddenValleyModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiddenValleyModel & operator=(const HiddenValleyModel &);

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
   *  Symmetry group for the new unbroken symmetry
   */
  //@{
  /**
   *  Type of group
   */
  GroupType groupType_;

  /**
   *  Order of the group
   */
  unsigned int Nc_;
  //@}

  /**
   *  Number of fermions
   */
  unsigned int Nf_;

  /**
   *  Charges etc under the new \f$U(1)\f$ group
   */
  //@{
  /**
   *  Charge of the left-handed quarks
   */
  double qL_;

  /**
   *  Charge of the right-handed up-type quarks
   */
  double uR_;

  /**
   *  Charge of the right-handed down-type quarks
   */
  double dR_;

  /**
   *  Charge of the left-handed leptons
   */
  double lL_;

  /**
   *  Charge of the right-handed charged leptons
   */
  double lR_;

  /**
   *  Coupling
   */
  double gPrime_;

  /**
   *  the vertex
   */
  AbstractFFVVertexPtr FFZPVertex_;

  /**
   *  The couplings of the quirks
   */
  vector<double> qCharge_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiddenValleyModel. */
template <>
struct BaseClassTrait<Herwig::HiddenValleyModel,1> {
  /** Typedef of the first base class of HiddenValleyModel. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiddenValleyModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiddenValleyModel>
  : public ClassTraitsBase<Herwig::HiddenValleyModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiddenValleyModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiddenValleyModel is implemented. It may also include several, space-separated,
   * libraries if the class HiddenValleyModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwHiddenValleyModel.so"; }
};

ThePEG_DECLARE_POINTERS(Herwig::HiddenValleyModel,HiddenValleyPtr);

/** @endcond */

}

#endif /* HERWIG_HiddenValleyModel_H */
