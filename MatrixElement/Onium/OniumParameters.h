// -*- C++ -*-
#ifndef Herwig_OniumParameters_H
#define Herwig_OniumParameters_H
//
// This is the declaration of the OniumParameters class.
//

#include "OniumParameters.fh"
#include "ThePEG/Interface/Interfaced.h"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

/**
 *  Enum to define the type of onium state
 */
enum OniumState {ccbar=0,bbbar=1,bcbar=2};
  
/**
 * The OniumParameters class stores the parameters for quarkonium production
 *
 * @see \ref OniumParametersInterfaces "The interfaces"
 * defined for OniumParameters.
 */
class OniumParameters: public Interfaced {

public:

  /**
   * The default constructor.
   */
  OniumParameters() : singletFromWaveFunction_(true),
		      R02_  (vector<vector<Energy3> >(3,vector<Energy3>())),
		      Rp02_ (vector<vector<Energy5> >(3,vector<Energy5>())),
		      Rpp02_(vector<vector<Energy7> >(3,vector<Energy7>())),
		      O1_S_(vector<vector<vector<Energy3> > >(3,vector<vector<Energy3> >())),
		      O1_P_(vector<vector<vector<Energy5> > >(3,vector<vector<Energy5> >())),
		      O1_D_(vector<vector<vector<Energy7> > >(3,vector<vector<Energy7> >()))
  {}

public:

  // Get the singlet matrix element
  template <unsigned int L>
  ThePEG::Qty<std::ratio<0,1>, std::ratio<3+2*L,1>, std::ratio<0,1>> inline
  singletME(OniumState type,unsigned int n, unsigned int S, unsigned int J);

public:

  /**
   *   Set the values of the wavefunction at the origin
   */
  string setWaveFunction(string arg);

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
  OniumParameters & operator=(const OniumParameters &) = delete;

private :

  /**
   *  Calculate the singlet matrix elements from the wavefunctions
   */
  bool singletFromWaveFunction_;

  /**
   *  Wavefunctions
   */
  //@{
  /**
   * \f$|R(0)|^2\f$ for the \f$s\f$-wave states
   */
  vector<vector<Energy3> > R02_;
  
  /**
   * \f$|R'(0)|^2\f$ for the \f$p\f$-wave states
   */
  vector<vector<Energy5> > Rp02_;
  
  /**
   * \f$|R''(0)|^2\f$ for the \f$d\f$-wave states
   */
  vector<vector<Energy7> > Rpp02_;
  //@}
  
  /**
   *  Singlet matrix elements
   */
  //@{
  /**
   * \f$|R(0)|^2\f$ for the \f$s\f$-wave states
   */
  vector<vector<vector<Energy3> > > O1_S_;
  
  /**
   * \f$|R'(0)|^2\f$ for the \f$p\f$-wave states
   */
  vector<vector<vector<Energy5 > > > O1_P_;
  
  /**
   * \f$|R''(0)|^2\f$ for the \f$d\f$-wave states
   */
  vector<vector<vector<Energy7> > > O1_D_;
  //@}
};

// Get the singlet matrix element
// s-wave
template <>
ThePEG::Qty<std::ratio<0,1>, std::ratio<3,1>, std::ratio<0,1>>
inline OniumParameters::singletME<0>(OniumState type, unsigned int n, unsigned int S, unsigned int J) {
  assert(O1_S_[type].size()>=n);
  assert(S==J && S<=1);
  return O1_S_[type][n-1][J];
}
// p-wave
template <>
ThePEG::Qty<std::ratio<0,1>, std::ratio<5,1>, std::ratio<0,1>>
inline OniumParameters::singletME<1>(OniumState type, unsigned int n, unsigned int S, unsigned int J) {
  assert(O1_P_[type].size()>=n);
  assert(S<=1&&J<=2);
  if(S==0) {
    assert(J==1);
    return O1_P_[type][n-1][0];
  }
  else {
    return O1_P_[type][n-1][J+1];
  }
}
// d-wave
template <>
ThePEG::Qty<std::ratio<0,1>, std::ratio<7,1>, std::ratio<0,1>>
inline OniumParameters::singletME<2>(OniumState type, unsigned int n, unsigned int S, unsigned int J) {
  assert(O1_D_[type].size()>=n);
  assert(S<=1&&J>0&&J<=3);
  if(S==0) {
    assert(J==2);
    return O1_D_[type][n-1][0];
  }
  else {
    return O1_D_[type][n-1][J];
  }
}
}

#endif /* Herwig_OniumParameters_H */
