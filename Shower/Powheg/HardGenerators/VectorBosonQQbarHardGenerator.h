// -*- C++ -*-
#ifndef HERWIG_VectorBosonQQbarHardGenerator_H
#define HERWIG_VectorBosonQQbarHardGenerator_H
//
// This is the declaration of the VectorBosonQQbarHardGenerator class.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "Herwig++/Shower/Powheg/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "VectorBosonQQbarHardGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/LorentzRotation.h"



namespace Herwig {

using namespace ThePEG;
using namespace std;
/**
 * Here is the documentation of the VectorBosonQQbarHardGenerator class.
 *
 * @see \ref VectorBosonQQbarHardGeneratorInterfaces "The interfaces"
 * defined for VectorBosonQQbarHardGenerator.
 */
class VectorBosonQQbarHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  inline VectorBosonQQbarHardGenerator();

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual NasonTreePtr generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  bool getEvent();
 
  /**
   * Constructs the post-emission momenta of q, qbar, g
   */
  void constructVectors();

  /** Returns the value of the radiative cross section (R(v,r)/B(v))
   *for the current (_xq,_xqb)
   */
  double getResult();

  /**
   *Returns true if we are within the allowed phase space of emission
   */   
  inline bool inRange();
 

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VectorBosonQQbarHardGenerator> 
  initVectorBosonQQbarHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VectorBosonQQbarHardGenerator & operator=(const VectorBosonQQbarHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  The alphaS for the overestimate of the true distribution
   */
  double _alphaS_max;  

  /**
   *  The dalitz variables (xq,xqb,xg=2-xq-xqb). These are
   *  the COM energies of the q,qb,g divided by 0.5*sqrt(_s).
   */
  double _xq;
  double _xqb;
  double _xg;

  /**
   *  The COM 3-momenta and of the q,qb,g divided by 0.5*sqrt(_s).
   *  These are used in checking phase space constraints as well
   *  as constructing emissions. This is beta (velocity) times 
   *  xq/xqb/xg hence the variable name.
   */
  double _b_xq;
  double _b_xqb;
  double _b_xg;

  /**
   *  _rt_mlambda/(2.*_b_xi*_b_xj) gives the sine of the angle 
   *  between i & j in the COM frame. It is equal to the square
   *  root of the Kallen function taking as its 3 arguments
   *  the rescaled 3-momenta _b_xq^2, _b_xqb^2, _b_xg^2. It 
   *  it is therefore importantly used in calculating relative pT's.
   */
  double _rt_mlambda;

  /**
   *  The gluon mass (Q_g) i.e. the cut off on pt, assuming massless quarks.
   */
  Energy _Qg;

  //  radiative variables (pt,y)
  double _y;
  Energy _pt;

  //com energy
  Energy2 _s;
 
  // The phi angle of the rotation of emitter-gluon plane 
  double _phi;

  // iemit = 0 quark emission: =1 antiquark emission
  int _iemitter;
  int _ispectator;

  HistogramPtr _hyb;
  HistogramPtr _hplow;
  HistogramPtr _hphigh;
  HistogramPtr _hy;
  HistogramPtr _hthrust;
  HistogramPtr _hthrustlow;
  HistogramPtr _hmass;
  HistogramPtr _hptVeto;
  HistogramPtr _hptVetoLow;
  /**
   *  vector of points for scatter plot
   */
  std::vector<double> _ptplot;
  std::vector<double> _yplot;
  std::vector<double> _xqplot;
  std::vector<double> _xqbplot;

  // The quark momenta and data pointers
  vector<Lorentz5Momentum> _quark;
  vector<tcPDPtr> _partons;
  PPtr _boson;
  // The gluon momentum and data pointer
  Lorentz5Momentum _g;
  tcPDPtr _gluon_data;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VectorBosonQQbarHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VectorBosonQQbarHardGenerator,1> {
  /** Typedef of the first base class of VectorBosonQQbarHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VectorBosonQQbarHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VectorBosonQQbarHardGenerator>
  : public ClassTraitsBase<Herwig::VectorBosonQQbarHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VectorBosonQQbarHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VectorBosonQQbarHardGenerator is implemented. 
   * It may also include several, space-separated,
   * libraries if the class VectorBosonQQbarHardGenerator depends
   * on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "VectorBosonQQbarHardGenerator.icc"

#endif /* HERWIG_VectorBosonQQbarHardGenerator_H */
