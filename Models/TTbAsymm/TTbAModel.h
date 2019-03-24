// -*- C++ -*-
#ifndef HERWIG_TTbAModel_H
#define HERWIG_TTbAModel_H
//
// This is the declaration of the TTbAModel class.
//

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "TTbAModel.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Here is the documentation of the TTbAModel class.
 *
 * @see \ref TTbAModelInterfaces "The interfaces"
 * defined for TTbAModel.
 */
class TTbAModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  TTbAModel();
 
  /** @name Vertices */
  //@{
 /**
   * Pointer to the object handling W prime vertex.
   */
  tAbstractFFVVertexPtr vertexWPTD() const {return _theWPTDVertex;}

  /**
   * Pointer to the object handling Z prime vertex.
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
   * Return the W prime top-down left-handed coupling
   */
  double _cWPTD_left() const {return _gWPTD_L;}

  /**
   * Return the W prime top-down right-handed coupling
   */
  double _cWPTD_right() const {return _gWPTD_R;}

  /**
   * Return the Z prime top-up left-handed coupling
   */
  double _cZPTU_left() const {return _gZPTU_L;}

  /**
   * Return the Z prime top-up right-handed coupling
   */
  double _cZPTU_right() const {return _gZPTU_R;}

  /**
   * Return the Z prime up-upbar left-handed coupling
   */
  double _cZPUU_left() const {return _gZPUU_L;}

  /**
   * Return the Z prime up-upbar right-handed coupling
   */
  double _cZPUU_right() const {return _gZPUU_R;}

    /**
   * Return the Z prime charm-charmbar left-handed coupling
   */
  double _cZPCC_left() const {return _gZPCC_L;}

  /**
   * Return the Z prime charm-charmbar right-handed coupling
   */
  double _cZPCC_right() const {return _gZPCC_R;}

  /**
   * Return the axigluon q-qbar left-handed coupling
   */
  double _cAGQQ_left() const {return _gAGQQ_L;}

  /**
   * Return the axigluon q-qbar right-handed coupling
   */
  double _cAGQQ_right() const {return _gAGQQ_R;}


  /**
   * Return the axigluon t-tbar left-handed coupling
   */
  double _cAGTT_left() const {return _gAGTT_L;}

  /**
   * Return the axigluon t-tbar right-handed coupling
   */
  double _cAGTT_right() const {return _gAGTT_R;}

  /**
   * Return the alphaX value of the SU(2)_X model
   */
  double _alphaX_value() const {return _alphaXparam;}

  /**
   * Return the costheta misalignment value of the SU(2)_X model
   */
  double _costhetaX_value() const {return _costhetaXparam;}

  /**
   * Return the selected model id
   */
  int _model() const {return _modelselect;}



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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TTbAModel & operator=(const TTbAModel &) = delete;

  
  /**
   * Pointer to the object handling the Wp to Top Down vertex.
   */
  AbstractFFVVertexPtr  _theWPTDVertex;

  
  /**
   * Pointer to the object handling the Zp to Quark-antiQuark vertex.
   */
  AbstractFFVVertexPtr  _theZPQQVertex;

  /**
   * Pointer to the object handling the Ag to Quark-antiQuark vertex.
   */
  AbstractFFVVertexPtr  _theAGQQVertex;

  /**
   * Pointer to the object handling the SU(2)_X vertex.
   */
  AbstractFFVVertexPtr  _theSU2XVertex;




 /**
   *  W prime coupling to top-down (left-handed)
   */
  double _gWPTD_L;
  

  /**
   *  W prime coupling to top-down (right-handed)
   */
  double _gWPTD_R;


   /**
   *  Z prime coupling to top-up (left-handed)
   */
  double _gZPTU_L;
  

  /**
   *  Z prime coupling to top-up (right-handed)
   */
  double _gZPTU_R;


   /**
   *  Z prime coupling to up-upbar (left-handed)
   */
  double _gZPUU_L;
  

  /**
   *  Z prime coupling to up-upbar (right-handed)
   */
  double _gZPUU_R;


   /**
   *  Z prime coupling to charm-charmbar (left-handed)
   */
  double _gZPCC_L;
  

  /**
   *  Z prime coupling to charm-charmbar (right-handed)
   */
  double _gZPCC_R;



   /**
   *  Axigluon coupling to q-qbar (left-handed)
   */
  double _gAGQQ_L;
  

  /**
   *  Axigluon coupling to q-qbar (right-handed)
   */
  double _gAGQQ_R;


  /**
   *  Axigluon coupling to t-tbar (left-handed)
   */
  double _gAGTT_L;
  

  /**
   *  Axigluon coupling to t-tbar (right-handed)
   */
  double _gAGTT_R;

 /**
   *  SU(2)_X alpha_X parameter
   */
  double _alphaXparam;
  
 /**
   *  SU(2)_X costheta misalignment angle
   */
  double _costhetaXparam;
  
 
  /** 
   * Model selector
   */
  int _modelselect;


};

}

#endif /* HERWIG_TTbAModel_H */
