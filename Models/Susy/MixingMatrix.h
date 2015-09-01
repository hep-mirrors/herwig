// -*- C++ -*-
//
// MixingMatrix.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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

  /**
   * A vector of mixing elements.
   */
  typedef vector<MixingElement> MixingVector;

  /**
   * The size of the matrix 
   */
  typedef pair<unsigned int, unsigned int> MatrixSize;
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

  /** Exception class to indicate problem with mixing matrix .*/
  class MixingMatrixError : public Exception {};

public:

  /** @name Constructors */
  //@{
  /**
   * Constructor that takes a mixing matrix and a vector id's as arguments 
   * @param mix Mixing matrix
   * @param ids The ids of the mixing sparticles
   */
  MixingMatrix(const CMatrix & mix,const vector<long> & ids) : 
    mixingMatrix_(mix),ids_(ids), size_(make_pair(mix.size(),mix[0].size())) 
  {}

  /**
   * Contructor that initializes size of matrix
   */
  MixingMatrix(unsigned int row, unsigned int col) :
    mixingMatrix_(row,vector<Complex>(col,Complex(0.,0.))), size_(row,col)
  {}

  /**
   * Standard Constructor.
   */
  MixingMatrix() {}
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
  void setMatrix(const CMatrix & mixing)  {
    mixingMatrix_ = mixing;
    size_ = make_pair(mixing.size(),mixing[0].size());
  }

  /**
   *Get the mixing matrix
   */
  CMatrix getMatrix() const {return mixingMatrix_;}

  /**
   * Set the vector containing mixing particles codes
   * @param mixingCodes vector containing PDG codes for mixing particles
   */
  void setIds(const vector<long> & mixingCodes)  {
    if(mixingCodes.size() != size_.first) {
      ostringstream codes;
      for(unsigned int ix=0;ix<mixingCodes.size();++ix)
	codes << mixingCodes[ix] << " ";
      throw MixingMatrixError() << "MixingMatrix::setIds() - The number "
				<< "of PDG codes (" << mixingCodes.size()
				<< ") does not match the size of the "
				<< "matrix (" << size_.first
				<< ")"
				<< "Ids are " << codes.str() 
				<< Exception::warning;
      return;
    }
    ids_ = mixingCodes;
  }
  
  /**
   * Get the vector containing mixing particles codes
   */
  const vector<long> & getIds() const {return ids_;}
  //@}
  
  /**
   * Multiply row corresponding to id by \f$i\f$
   * @param id PDG code of particle
   */
  void adjustPhase(long id);
    
  /**
   * Access element of matrix
   */  
  const Complex operator()(unsigned int row, unsigned int col) const {
    return mixingMatrix_.at(row).at(col);
  }

  /**
   * Set element of matrix
   */
  Complex & operator()(unsigned int row, unsigned int col) {
    return mixingMatrix_.at(row).at(col);
  }
  
  /**
   * Add a PDG code to the stored vector
   */
  void addCode(long id) {
    if(ids_.size() >= size_.first) {
      throw MixingMatrixError() << "MixingMatrix::addCode() - Trying to add a"
				<< "PDG code but the vector already contains the "
				<< "same number as the matrix size " 
				<< Exception::warning;
      return;
    }
    ids_.push_back(id);
  }
  
  /**
   * Return the size of the mixing matrix
   */
  MatrixSize size() const {return size_;}
  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MixingMatrix & operator=(const MixingMatrix &);

  /**
   * The mixing matrix
   */
  CMatrix mixingMatrix_;

  /**
   * The PDG codes of the mixing particles
   */
  vector<long> ids_;

  /**
   * Size of matrix
   */
  pair<unsigned int,unsigned int> size_;
  
  /**
   * Print the matrix to the stream
   */
  friend ostream & operator<<(ostream & os,const MixingMatrix & mix);
};

  /**
   * Output operator for the MixingMatrix
   */
  ostream & operator<<(ostream &,const MixingMatrix &);
}

#endif /* HERWIG_MixingMatrix_H */
