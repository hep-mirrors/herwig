// -*- C++ -*-
#ifndef Herwig_RPVSSSVertex_H
#define Herwig_RPVSSSVertex_H
//
// This is the declaration of the RPVSSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/SSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the RPVSSSVertex class.
 *
 * @see \ref RPVSSSVertexInterfaces "The interfaces"
 * defined for RPVSSSVertex.
 */
class RPVSSSVertex: public Helicity::SSSVertex {

public:

  /**
   * The default constructor.
   */
  RPVSSSVertex();

  /**
   * Calculate the coupling for the vertex
   * @param q2 The scale to at which evaluate the coupling.
   * @param particle1 The first particle in the vertex.
   * @param particle2 The second particle in the vertex.
   * @param particle3 The third particle in the vertex.
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr particle2,
			   tcPDPtr particle3);

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
  RPVSSSVertex & operator=(const RPVSSSVertex &);

private:

  /**
   *  Which types of interaction to include
   */
  unsigned int interactions_;

  /**
   * The scale at which the coupling was last evaluated.  
   */
  Energy2 q2Last_;




  //////////////   FROM HSFSF ///////////////////////////////////////////////////

//   /** @name Functions to calculate the coupling based on the sfermion type. */
//   //@{
//   /**
//    * Calculate the coupling for the first higgs 
//    * @param higgs The ID of the higgs
//    * @param smID The ID of the SM particle to which it is a partner.
//    * @param alpha The mass eigenstate of an sfermion
//    * @param beta  The mass eigenstate of the other sfermion
//    */
//   void downSF(long higgs, long smID, unsigned int alpha, unsigned int beta);
  
//   /**
//    * Calculate the coupling for the second higgs 
//    * @param higgs The ID of the higgs
//    * @param smID The ID of the SM particle to which it is a partner.
//    * @param alpha The mass eigenstate of an sfermion
//    * @param beta  The mass eigenstate of the other sfermion
//    */
//   void upSF(long higgs, long smID, unsigned int alpha, unsigned int beta);
  
//   /**
//    * Calculate the coupling for the third higgs 
//    * @param higgs The ID of the higgs
//    * @param smID The ID of the SM particle to which it is a partner.
//    * @param alpha The mass eigenstate of an sfermion
//    * @param beta  The mass eigenstate of the other sfermion
//    */
//   void leptonSF(long higgs, long smID, unsigned int alpha, unsigned int beta);
  
//   /**
//    *  Calculate the coupling for the charged higgs 
//    * @param id1 The ID of the first sfermion
//    * @param id2 The ID of the second sfermion
//    */
//   void chargedHiggs(long id1, long id2);
  
//   //@}
// private:
  
//   /**
//    * A vector containing pointers to the mixing matrices, 0 = stop, 
//    * 1 = sbottom, 2 = stau 
//    */
//   MMPVector theMix;

//   /**
//    * A vector containing the trilinear couplings, quarks then leptons
//    */
//   vector<complex<Energy> > theTriC;
    
//   /**
//    * The value of \f$\sin\alpha\f$.
//    */
//   double theSinA;

//   /**
//    * The value of \f$\cos\alpha\f$.
//    */
//   double theCosA;

//   /**
//    * The value of \f$\sin\beta\f$.
//    */
//   double theSinB;

//   /**
//    * The value of \f$\cos\beta\f$.
//    */
//   double theCosB;

//   /**
//    * The value of \f$\tan\beta\f$. 
//    */
//   double theTanB;

//   /**
//    * The value of \f$\sin(\alpha + \beta)\f$.
//    */
//   double theSinAB;

//   /**
//    * The value of \f$\cos(\alpha + \beta)\f$.
//    */
//   double theCosAB;

//   /**
//    * The mass of the \f$W\f$. 
//    */
//   Energy theMw;

//   /**
//    * The mass of the \f$Z\f$. 
//    */
//   Energy theMz;
  
//   /**
//    * The \f$\mu\f$ parameter. 
//    */
//   Energy theMu;

//   /**
//    * The value of \f$\sin\theta_W\f$
//    */
//   double theSw;

//   /**
//    * The value of \f$\cos\theta_W\f$
//    */
//   double theCw;

//   /**
//    * The value of the coupling when it was last evaluated
//    */
//   complex<Energy> theCoupLast;
  
//   /**
//    * The value of g coupling when it was last evaluated
//    */
//   double thegLast;
  
//   /**
//    * The ID of the higgs when the vertex was last evaluated 
//    */
//   long theHLast;
  
//   /**
//    * The ID of the first sfermion when the vertex was last evaluated 
//    */
//   long theSF1Last;
  
//   /**
//    * The ID of the second sfermion when the vertex was last evaluated 
//    */
//   long theSF2Last;
  
////////////////////////// FROM HHH ////////////////////////////////////////////

// private:
  
//   /**
//    * The mass of the \f$W\f$.
//    */
//   Energy theMw;

//   /**
//    * The factor \f$ \frac{m_Z}{\sin2\theta_W} \f$
//    */
//   Energy theZfact;

//   /**
//    * The value of \f$\sin\theta_W\f$
//    */
//   double theSw;

//   /**
//    * The value of \f$ \sin(\beta + \alpha) \f$.
//    */
//   double theSbpa;

//   /**
//    * The value of \f$ \cos(\beta + \alpha) \f$.
//    */
//   double theCbpa;

//   /**
//    * The value of \f$ \sin(\beta - \alpha) \f$.
//    */
//   double theSbma;

//   /**
//    * The value of \f$ \cos(\beta - \alpha) \f$.
//    */
//   double theCbma;
  
//   /**
//    * The value of \f$ \sin 2\alpha \f$.
//    */
//   double theS2a;

//   /**
//    * The value of \f$ \cos 2\alpha \f$.
//    */
//   double theC2a;

//   /**
//    * The value of \f$ \sin 2\beta \f$.
//    */
//   double theS2b;

//   /**
//    * The value of \f$ \cos 2\beta \f$.
//    */
//   double theC2b;
  
//   /**
//    * The value of \f$ \sqrt{4\pi\alpha}\f$ when it was last evaluated.
//    */
//   double theElast;

};

}

#endif /* Herwig_RPVSSSVertex_H */
