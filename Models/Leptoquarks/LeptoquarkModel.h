// -*- C++ -*-
#ifndef HERWIG_LeptoquarkModel_H
#define HERWIG_LeptoquarkModel_H
//
// This is the declaration of the LeptoquarkModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
#include "LeptoquarkModel.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the LeptoquarkModel class.
 *
 * @see \ref LeptoquarkModelInterfaces "The interfaces"
 * defined for LeptoquarkModel.
 */
class LeptoquarkModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  LeptoquarkModel();
 
  /** @name Vertices */
  //@{
  /**
   * Pointer to the object handling S0S0barg vertex.
   */
  tAbstractVSSVertexPtr   vertexSLQSLQG() const {return _theSLQSLQGVertex;}
  
  /**
   * Pointer to the object handling the S0S0bargg vertex.
   */
  tAbstractVVSSVertexPtr   vertexSLQSLQGG() const {return _theSLQSLQGGVertex;}

  /**
   * Pointer to the object handling the S0ql vertex.
   */
  tAbstractFFSVertexPtr   vertexSLQFF() const {return _theSLQFFVertex;}

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

 /**
   * Return the overall fermion coupling
   */
  double cfermion() const {return _CouplFF;}

  /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (S0)
   */
  double cleft() const {return _leftcoup;}

 /**
   * Return the coupling of the leptoquark to right-handed leptons + left-handed quarks (S0)
   */
  double cright() const {return _rightcoup;}

  /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (S1 triplet)
   */
  double cleft1() const {return _leftcoup1;}

  /**
   * Return the coupling of the leptoquark to right-handed leptons
   * + left-handed quarks (~S0)
   */
  double crighttilde() const {return _rightcouptilde;}

 /**
  * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (S1/2)
   */  
  double cleft12() const {return _leftcoup12;}

  /**
   * Return the coupling of the leptoquark to right-handed leptons + left-handed quarks (S1/2)
   */
  double cright12() const {return _rightcoup12;}

  /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (S1/2)
   */
  double cleft12tilde() const {return _leftcoup12t;}

 /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (S0)
   */
  double dcleft() const {return _dleftcoup;}
 /**
   * Return the coupling of the leptoquark to right-handed leptons + left-handed quarks (dS0)
   */
  double dcright() const {return _drightcoup;}

  /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (dS1 triplet)
   */

  double dcleft1() const {return _dleftcoup1;}

  /**
   * Return the coupling of the leptoquark to right-handed leptons
   * + left-handed quarks (~dS0)
   */
  double dcrighttilde() const {return _drightcouptilde;}

 /**
  * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (dS1/2)
   */  
  double dcleft12() const {return _dleftcoup12;}

  /**
   * Return the coupling of the leptoquark to right-handed leptons + left-handed quarks (dS1/2)
   */
  double dcright12() const {return _drightcoup12;}

  /**
   * Return the coupling of the leptoquark to left-handed leptons + right-handed quarks (dS1/2)
   */
  
  double dcleft12tilde() const {return _dleftcoup12t;}
  /**
   * Suppression scale for derivatively coupled scalar leptoquarks
   */
  
  Energy fscale() const {return _derivscalef;}





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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<LeptoquarkModel> initLeptoquarkModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LeptoquarkModel & operator=(const LeptoquarkModel &) = delete;

  
  /**
   * Pointer to the object handling the G to SLQ SLQ vertex.
   */
  AbstractVSSVertexPtr  _theSLQSLQGVertex;

  /**
   * Pointer to the object handling the GG to SLQ SLQ vertex.
   */
  AbstractVVSSVertexPtr  _theSLQSLQGGVertex;


  /**
   * Pointer to the object handling the SLQ to FF vertex.
   */
  AbstractFFSVertexPtr  _theSLQFFVertex;



  /**
   *  Overall coupling to fermions
   */
  double _CouplFF;


  /**
   *  Overall coupling to left-handed leptons (S0)
   */
  double _leftcoup;

  /**
   *  Overall coupling to right-handed leptons (S0)
   */
  double _rightcoup;

  /**
   *  Overall coupling to left-handed leptons (~S0)
   */
  double _rightcouptilde;
  
  /**
   *  Overall coupling to left-handed leptons (S1 triplet)
   */
  double _leftcoup1;

  /**
   *  Overall coupling to left-handed leptons (S1/2)
   */
  double _leftcoup12;

  /**
   *  Overall coupling to right-handed leptons (S1/2)
   */
  double _rightcoup12;

  /**
   *  Overall coupling to left-handed leptons (~S1/2)
   */
  double _leftcoup12t;

   /**
   *  Overall coupling to left-handed leptons (dS0)
   */
  double _dleftcoup;

  /**
   *  Overall coupling to right-handed leptons (dS0)
   */
  double _drightcoup;

  /**
   *  Overall coupling to left-handed leptons (~dS0)
   */
  double _drightcouptilde;
  
  /**
   *  Overall coupling to left-handed leptons (dS1 triplet)
   */
  double _dleftcoup1;

  /**
   *  Overall coupling to left-handed leptons (dS1/2)
   */
  double _dleftcoup12;

  /**
   *  Overall coupling to right-handed leptons (dS1/2)
   */
  double _drightcoup12;

  /**
   *  Overall coupling to left-handed leptons (~dS1/2)
   */
  double _dleftcoup12t;

  /**
   *  Suppression scale for derivatively coupled scalar leptoquarks, f
   */
  Energy _derivscalef;

  

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LeptoquarkModel. */
template <>
struct BaseClassTrait<Herwig::LeptoquarkModel,1> {
  /** Typedef of the first base class of LeptoquarkModel. */
  typedef Herwig::BSMModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LeptoquarkModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LeptoquarkModel>
  : public ClassTraitsBase<Herwig::LeptoquarkModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LeptoquarkModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LeptoquarkModel is implemented. It may also include several, space-separated,
   * libraries if the class LeptoquarkModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwLeptoquarkModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LeptoquarkModel_H */
