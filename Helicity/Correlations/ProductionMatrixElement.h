// -*- C++ -*-
#ifndef HERWIG_ProductionMatrixElement_H
#define HERWIG_ProductionMatrixElement_H
//
// This is the declaration of the ProductionMatrixElement class.

#include <ThePEG/Config/ThePEG.h>
#include <ThePEG/Utilities/ClassDescription.h>
#include <ThePEG/Helicity/RhoDMatrix.h>
// #include "ProductionMatrixElement.fh"
// #include "ProductionMatrixElement.xh"

namespace Herwig {
using ThePEG::Helicity::RhoDMatrix;
namespace Helicity {

using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The storage of the helicity amplitude expression for the matrix element 
 *  of a hard process. Two incoming particles and an arbitary number of 
 *  external particles are supported.
 *
 *  @see DecayMatrixElement
 *  @see RhoDMatrix
 *  @see HardVertex
 *
 *  \author Peter Richardson
 */

class ProductionMatrixElement: public Base {  
      
public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor for 2-1 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out \f$2S+1\f$ for the outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out);

  /**
   * Constructor for 2-2 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out1,int out2);

  /**
   * Constructor for 2-3 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out1,int out2,int out3);

  /**
   * Constructor for 2-4 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out1,int out2,int out3, int out4);

  /**
   * Constructor for 2-5 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   * @param out5 \f$2S+1\f$ for the fifth outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out1,int out2,int out3, int out4,
				 int out5);

  /**
   * Constructor for 2-6 scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out1 \f$2S+1\f$ for the first outgoing particle.
   * @param out2 \f$2S+1\f$ for the second outgoing particle.
   * @param out3 \f$2S+1\f$ for the third outgoing particle.
   * @param out4 \f$2S+1\f$ for the fourth outgoing particle.
   * @param out5 \f$2S+1\f$ for the fifth outgoing particle.
   * @param out6 \f$2S+1\f$ for the sixth outgoing particle.
   */
  inline ProductionMatrixElement(int in1,int in2,int out1,int out2,int out3, int out4,
				 int out5, int out6);

  /**
   * Constructor for 2-n scattering.
   * @param in1 \f$2S+1\f$ for the first incoming particle.
   * @param in2 \f$2S+1\f$ for the second incoming particle.
   * @param out A vector containing \f$2S+1\f$ for the outgoing particles.
   */
  inline ProductionMatrixElement(int in1,int in2,vector<int> out);

  /**
   * Default constructor.
   */
  inline ProductionMatrixElement();

  /**
   * Copy-constructor.
   */
  inline ProductionMatrixElement(const ProductionMatrixElement &);

  /**
   * Destructor.
   */
  virtual ~ProductionMatrixElement();
  //@}

public:
     
  /** @name Access to the spins of the particles. */
  //@{
  /**
   * Get the spins of the incoming particles particle
   * @return A vector containing \f$2S+1\f$ for the two incoming particles.
   */
  inline vector<int> inspin();

  /**
   * Get the spins of the outgoing particles.
   * @return A vector containing \f$2S+1\f$ for the outgoing particles.
   */
  inline vector<int> outspin();
  //@}  

public:

  /** @name Access to the individual helicity components. */
  //@{

  /**
   * Access the helicity components for a 2-1 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel The helicity of the outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel) const;

  /**
   * Access the helicity components for a 2-1 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel The helicity of the outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel);

  /**
   * Access the helicity components for a 2-2 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel1,int outhel2) const;

  /**
   * Access the helicity components for a 2-2 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel1,int outhel2);

  /**
   * Access the helicity components for a 2-3 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3) const;

  /**
   * Access the helicity components for a 2-3 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3);

  /**
   * Access the helicity components for a 2-4 scattering.  This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3,int outhel4) const;

  /**
   * Access the helicity components for a 2-4 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3, int outhel4);

  /**
   * Access the helicity components for a 2-5 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @param outhel5 The helicity of the fifth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3,int outhel4, int outhel5) const;

  /**
   * Access the helicity components for a 2-5 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @param outhel5 The helicity of the fifth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3, int outhel4, int outhel5);

  /**
   * Access the helicity components for a 2-6 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @param outhel5 The helicity of the fifth outgoing particle.
   * @param outhel6 The helicity of the sixth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3,int outhel4, int outhel5,
				int outhel6) const;

  /**
   * Access the helicity components for a 2-6 scattering. This method supplies
   * the component and allows it to be changed.
   * @param inhel1 The helicity of the first incoming particle.
   * @param inhel2 The helicity of the second incoming particle.
   * @param outhel1 The helicity of the first outgoing particle.
   * @param outhel2 The helicity of the second outgoing particle.
   * @param outhel3 The helicity of the third outgoing particle.
   * @param outhel4 The helicity of the fourth outgoing particle.
   * @param outhel5 The helicity of the fifth outgoing particle.
   * @param outhel6 The helicity of the sixth outgoing particle.
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (int inhel1,int inhel2,int outhel1,int outhel2,
				int outhel3, int outhel4, int outhel5, int outhel6);


  /**
   * Access the helicity components for a 2-6 scattering. This method supplies
   * the component but does not allow it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  inline Complex   operator () (vector<int> hel) const;

  /**
   * Access the helicity components for a 2-n scattering. This method supplies
   * the component and allows it to be changed.
   * @param hel The helicities of the incoming and outgoing particles
   * @return The matrix element for the given helicities.
   */
  inline Complex & operator () (vector<int> hel);
  //@}

public:

  /**
   * Calculate the decay matrix for an incoming particle.
   */
  RhoDMatrix calculateDMatrix(int,RhoDMatrix, vector<RhoDMatrix>);

  /**
   * Calculate the rho matrix for a given outgoing particle.
   */
  RhoDMatrix calculateRhoMatrix(int,RhoDMatrix,RhoDMatrix,
				vector<RhoDMatrix>);
  
public:

  /**
   * Reset the matrix element.
   */
  inline void reset(const ProductionMatrixElement &) const;

public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
private:
  
  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<ProductionMatrixElement> initProductionMatrixElement;
  
  /**
   * Private and non-existent assignment operator.
   */
  ProductionMatrixElement & operator=(const ProductionMatrixElement& );
  
private:
  
  /**
   * Set the size of the vector containing the matrix element.
   */
  inline void setMESize();
  
private:
  
  /**
   * Number of outgoing particles.
   */
  mutable int _nout;

  /**
   * Spin of the incoming particles as 2s+1.
   */
  mutable vector<int> _inspin;

  /**
   * Spins of the outgoing particles.
   */
  mutable vector<int> _outspin;

  /**
   * Storage of the matrix element, a vector is better for memory usage.
   */
  mutable vector<Complex> _matrixelement;

  /**
   * Constants needed to map the index of the vector to a helicity structure.
   */
  mutable vector<int> _constants;

};

}
}


namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of ProductionMatrixElement.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::ProductionMatrixElement,1> {
    /** Typedef of the base class of ProductionMatrixElement. */
    typedef Base NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::ProductionMatrixElement>
    : public ClassTraitsBase<Herwig::Helicity::ProductionMatrixElement> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::Helicity::ProductionMatrixElement"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwCorrelations.so"; }

  };
  
}

#include "ProductionMatrixElement.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ProductionMatrixElement.tcc"
#endif

#endif /* HERWIG_ProductionMatrixElement_H */
