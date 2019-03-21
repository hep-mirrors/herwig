// -*- C++ -*-
#ifndef HERWIG_ZprimeModel_H
#define HERWIG_ZprimeModel_H
//
// This is the declaration of the ZprimeModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ZprimeModel.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the ZprimeModel class.
 *
 * @see \ref ZprimeModelInterfaces "The interfaces"
 * defined for ZprimeModel.
 */
class ZprimeModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  ZprimeModel();
 
  /** @name Vertices */
  //@{
 
  /**
   * Pointer to the object handling Z prime quark-anti-quark vertex.
   */
  tAbstractFFVVertexPtr vertexZPQQ() const {return _theZPQQVertex;}

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
   * Return the Z prime top-up left-handed coupling
   */
  double _cZPTU_left() const {return _gZPTU_L;}

  /**
   * Return the Z prime top-up right-handed coupling
   */
  double _cZPTU_right() const {return _gZPTU_R;}

  /**
   * Return the Z prime d-dbar left-handed coupling
   */
  double _cZPDD_left() const {return _gZPDD_L;}

  /**
   * Return the Z prime d-dbar right-handed coupling
   */
  double _cZPDD_right() const {return _gZPDD_R;}

  /**
   * Return the Z prime top-anti-top left-handed coupling
   */
  double _cZPTT_left() const {return _gZPTT_L;}

  /**
   * Return the Z prime top-anti-top right-handed coupling
   */
  double _cZPTT_right() const {return _gZPTT_R;}

  /**
   * Return the Z prime u-ubar left-handed coupling
   */
  double _cZPUU_left() const {return _gZPUU_L;}

  /**
   * Return the Z prime u-ubar right-handed coupling
   */
  double _cZPUU_right() const {return _gZPUU_R;}

    /**
   * Return the Z prime c-cbar left-handed coupling
   */
  double _cZPCC_left() const {return _gZPCC_L;}

  /**
   * Return the Z prime c-cbar right-handed coupling
   */
  double _cZPCC_right() const {return _gZPCC_R;}

  /**
   * Return the Z prime b-bbar left-handed coupling
   */
  double _cZPBB_left() const {return _gZPBB_L;}

  /**
   * Return the Z prime b-bbar right-handed coupling
   */
  double _cZPBB_right() const {return _gZPBB_R;}
  
    /**
   * Return the Z prime s-sbar left-handed coupling
   */
  double _cZPSS_left() const {return _gZPSS_L;}

  /**
   * Return the Z prime c-cbar right-handed coupling
   */
  double _cZPSS_right() const {return _gZPSS_R;}


  /**
   * Return the Z prime e+e- left-handed coupling
   */
  double _cZPee_left() const {return _gZPee_L;}

  /**
   * Return the Z prime e+e- right-handed coupling
   */
  double _cZPee_right() const {return _gZPee_R;}

  /**
   * Return the Z prime mu+mu- left-handed coupling
   */
  double _cZPmm_left() const {return _gZPmm_L;}

  /**
   * Return the Z prime mu+mu- right-handed coupling
   */
  double _cZPmm_right() const {return _gZPmm_R;}
  
    /**
   * Return the Z prime tau+tau- left-handed coupling
   */
  double _cZPtt_left() const {return _gZPtt_L;}

  /**
   * Return the Z prime tau+tau- right-handed coupling
   */
  double _cZPtt_right() const {return _gZPtt_R;}

  /**
   * Return the Z prime nu_e nu_ebar left-handed coupling
   */
  double _cZPnuenue_left() const {return _gZPnuenue_L;}

  /**
   * Return the Z prime nu_e nu_ebar right-handed coupling
   */
  double _cZPnuenue_right() const {return _gZPnuenue_R;}

  /**
   * Return the Z prime nu_mu nu_mubar left-handed coupling
   */
  double _cZPnumnum_left() const {return _gZPnumnum_L;}

  /**
   * Return the Z prime nu_mu nu_mubar right-handed coupling
   */
  double _cZPnumnum_right() const {return _gZPnumnum_R;}
  
    /**
   * Return the Z prime nu_tau nu_taubar left-handed coupling
   */
  double _cZPnutnut_left() const {return _gZPnutnut_L;}

  /**
   * Return the Z prime nu_tau nu_taubar right-handed coupling
   */
  double _cZPnutnut_right() const {return _gZPnutnut_R;}



  /**
   * Return the overall coupling of the Z prime to quark-anti-quark
   */
  double _cZPoverallCoup() const {return _ZPoverall;}



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
  static ClassDescription<ZprimeModel> initZprimeModel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ZprimeModel & operator=(const ZprimeModel &) = delete;

  
  /**
   * Pointer to the object handling the Zp to Quark-antiQuark vertex.
   */
  AbstractFFVVertexPtr  _theZPQQVertex;

    /**
   *  Z prime coupling to u-ubar (left-handed)
   */
  double _gZPUU_L;
  

  /**
   *  Z prime coupling to u-ubar (right-handed)
   */
  double _gZPUU_R;
 /**
   *  Z prime coupling to d-dbar (left-handed)
   */
  double _gZPDD_L;
  

  /**
   *  Z prime coupling to d-dbar (right-handed)
   */
  double _gZPDD_R;

 /**
   *  Z prime coupling to c-cbar (left-handed)
   */
  double _gZPCC_L;
  

  /**
   *  Z prime coupling to c-cbar (right-handed)
   */
  double _gZPCC_R;

 
 /**
   *  Z prime coupling to s-sbar (left-handed)
   */
  double _gZPSS_L;
  

  /**
   *  Z prime coupling to s-sbar (right-handed)
   */
  double _gZPSS_R;



  
 
/**
   *  Z prime coupling to b-bbar (left-handed)
   */
  double _gZPBB_L;
  

  /**
   *  Z prime coupling to b-bbar (right-handed)
   */
  double _gZPBB_R;

 /**
   *  Z prime coupling to top-up (left-handed)
   */
  double _gZPTU_L;
  

  /**
   *  Z prime coupling to top-up (right-handed)
   */
  double _gZPTU_R;

 /**
   *  Z prime coupling to top-anti-top (left-handed)
   */
  double _gZPTT_L;
  

  /**
   *  Z prime coupling to top-anti-top (right-handed)
   */
  double _gZPTT_R;



   /**
   *  Z prime coupling to e+e- (left-handed)
   */
  double _gZPee_L;
  

  /**
   *  Z prime coupling to e+e- (right-handed)
   */
  double _gZPee_R;

 
/**
   *  Z prime coupling to mu+mu- (left-handed)
   */
  double _gZPmm_L;
  

  /**
   *  Z prime coupling to mu+mu- (right-handed)
   */
  double _gZPmm_R;

 /**
   *  Z prime coupling to tau+tau- (left-handed)
   */
  double _gZPtt_L;
  

  /**
   *  Z prime coupling to tau+tau- (right-handed)
   */
  double _gZPtt_R;


   /**
   *  Z prime coupling to nu_e nu_ebar (left-handed)
   */
  double _gZPnuenue_L;
  

  /**
   *  Z prime coupling to nu_e nu_ebar (right-handed)
   */
  double _gZPnuenue_R;

 
/**
   *  Z prime coupling to nu_mu nu_mubar (left-handed)
   */
  double _gZPnumnum_L;
  

  /**
   *  Z prime coupling to nu_mu nu_mubar (right-handed)
   */
  double _gZPnumnum_R;

 /**
   *  Z prime coupling to nu_tau nu_taubar (left-handed)
   */
  double _gZPnutnut_L;
  

  /**
   *  Z prime coupling to nu_tau nu_taubar (right-handed)
   */
  double _gZPnutnut_R;


 /**
   *  Z prime overall coupling
   */
  double _ZPoverall;
  



};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ZprimeModel. */
template <>
struct BaseClassTrait<Herwig::ZprimeModel,1> {
  /** Typedef of the first base class of ZprimeModel. */
  typedef Herwig::StandardModel NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ZprimeModel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ZprimeModel>
  : public ClassTraitsBase<Herwig::ZprimeModel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ZprimeModel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ZprimeModel is implemented. It may also include several, space-separated,
   * libraries if the class ZprimeModel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwZprimeModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ZprimeModel_H */
