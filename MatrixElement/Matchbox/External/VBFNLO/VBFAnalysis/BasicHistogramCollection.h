// -*- C++ -*-
#ifndef ARNOLD_BasicHistogramCollection_H
#define ARNOLD_BasicHistogramCollection_H
//
// This is the declaration of the BasicHistogramCollection class.
//
#include "arHistogram.h"
#include <iostream>
#include "fastjet/PseudoJet.hh"
#include "ThePEG/Interface/Interfaced.h"
namespace arnold {

using namespace ThePEG;

class BasicHistogramCollection: public Interfaced{

public:
  //Constructor
  BasicHistogramCollection();

  //Destructor
  virtual ~BasicHistogramCollection();

  void writetofolder(const char * folder );

  void fill(const vector<fastjet::PseudoJet> & jets, const vector<fastjet::PseudoJet> & partons,const fastjet::PseudoJet &, const double w = 1.0);

  double A_phi_jets, A_phi_partons;

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

  int theNBins;

  bool doSmearing;
    
  arHistogramPtr histo[236];

  ar2DHistogramPtr histo2D[8];

  long correlarray[2][3];

  long Nj_delphi_25to75, Nj_delphi_0to25, Nj_delphi_75to100;
  long Np_delphi_25to75, Np_delphi_0to25, Np_delphi_75to100;

};

class SingleJetHistogramCollection {

public:
  //Constructor
  SingleJetHistogramCollection();

  //Destructor
  virtual ~SingleJetHistogramCollection();

  void writetofolder(const char * folder );

  void fill(const fastjet::PseudoJet &,const fastjet::PseudoJet & ,double w);



private:
  arHistogramPtr histo[5];
  

  
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMEBase. */
template <>
struct BaseClassTrait<arnold::BasicHistogramCollection,1> {
  /** Typedef of the first base class of BasicHistogramCollection. */
  typedef Interfaced NthBase;
}; 

/** This template specialization informs ThePEG about the name of
 *  the BasicHistogramCollection class and the shared object where it is defined. */
template <>
struct ClassTraits<arnold::BasicHistogramCollection>
  : public ClassTraitsBase<arnold::BasicHistogramCollection> {
  /** Return a platform-independent class name */
  static string className() { return "arnold::BasicHistogramCollection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FJAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class FJAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "VBFAnalysis.so"; }
};

/** @endcond */

}

#endif /* ARNOLD_BasicHistogramCollection_H */
