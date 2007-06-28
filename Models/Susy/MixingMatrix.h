// -*- C++ -*-
#ifndef HERWIG_MixingMatrix_H
#define HERWIG_MixingMatrix_H
//
// This is the declaration of the MixingMatrix class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Config/Complex.h"
#include "MixingMatrix.fh"
#include "ThePEG/Utilities/Exception.h"

namespace Herwig {
using namespace ThePEG;

/*@name Some convenient typedefs. */
//@{

  /**
   * A complex valued nested vector. 
   */
  typedef vector<vector<Complex> > CMatrix;
  
  /**
   * Struct for the elements of a mixing matrix
   */
  struct MixingElement{

    /**
     *  Constructor
     */
    MixingElement(unsigned int irow, unsigned int icol,Complex ivalue) 
      : row(irow), col(icol), value(ivalue) {}

    /**
     * row
     */
    unsigned int row;

    /**
     * column
     */
    unsigned int col;

    /**
     *  value
     */
    Complex value;
  };

//@}

/**
 * This class is desinged to store the mixing matrices needed for Susy
 * studies. The actual matrix is stored as a nested complex vector. It
 * also stores a vector of PDG codes correspoding to the mass states of 
 * mixing states.
 * 
 * @see Interfaced
 */
class MixingMatrix: public Interfaced {

public:

  /** @name Constructors */
  //@{
  /**
   * Constructor that takes a mixing matrix and a vector id's as arguments 
   * @param mix Mixing matrix
   * @param ids The ids of the mixing sparticles
   */
  inline MixingMatrix(const CMatrix & mix,const vector<long> & ids);

  /**
   * Contructor that initializes size of matrix
   */
  inline MixingMatrix(unsigned int col, unsigned int row);

  /**
   * Standard Constructor.
   */
   inline MixingMatrix();
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

  /** @name Get and Set functions. */
  //@{
  /**
   * Set the mixing matrix
   * @param mixing The Mixing matrix stored as nested complex vector
   */
  inline void setMatrix(const CMatrix & mixing);

  /**
   *Get the mixing matrix
   */
  inline CMatrix getMatrix() const;

  /**
   * Set the vector containing mixing particles codes
   * @param mixingCodes vector containing PDG codes for mixing particles
   */
  inline void setIds(const vector<long> & mixingCodes);
  
  /**
   * Get the vector containing mixing particles codes
   */
  inline const vector<long> & getIds() const;
  //@}
  
  /**
   * Multiply row corresponding to id by \f$i\f$
   * @param id PDG code of particle
   */
  void adjustPhase(long id);
    
  /**
   * Access element of matrix
   */  
  inline const Complex operator()(unsigned int row, 
				  unsigned int col) const;

  /**
   * Set element of matrix
   */
  inline Complex & operator()(unsigned int row,
			      unsigned int col) ;

  /**
   * Add a PDG code to the stored vector
   */
  inline void addCode(long id);
  
  /**
   * Return the size of the mixing matrix
   */
  inline pair<unsigned int,unsigned int> size() const;

  // /**
//    * Overloaded multiplication operator
//    */ 
//   inline MixingMatrix operator*(const MixingMatrix &) const;
  
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MixingMatrix> initMixingMatrix;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MixingMatrix & operator=(const MixingMatrix &);

  /**
   * The mixing matrix
   */
  CMatrix _theMixingMatrix;

  /**
   * The PDG codes of the mixing particles
   */
  vector<long> _theIds;

  /**
   * Size of matrix
   */
  pair<unsigned int,unsigned int> _theSize;
  
  /**
   * Print the matrix to the stream
   */
  friend ostream & operator<<(ostream & os,const MixingMatrix & mix);
};

  /**
   * Output operator for the MixingMatrix
   */
  ostream & operator<<(ostream &,const MixingMatrix &);

  /** Exception class to indicate problem with mixing matrix .*/
  class MixingMatrixError : public Exception {};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MixingMatrix. */
template <>
struct BaseClassTrait<Herwig::MixingMatrix,1> {
  /** Typedef of the first base class of MixingMatrix. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MixingMatrix class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MixingMatrix>
  : public ClassTraitsBase<Herwig::MixingMatrix> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MixingMatrix"; }
  /** Return the name of the shared library be loaded to get
   *  access to the MixingMatrix class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSusy.so"; }
};

/** @endcond */

}

#include "MixingMatrix.icc"

#endif /* HERWIG_MixingMatrix_H */
