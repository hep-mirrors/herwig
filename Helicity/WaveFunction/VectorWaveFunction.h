// -*- C++ -*-
#ifndef HERWIG_VectorWaveFunction_H
#define HERWIG_VectorWaveFunction_H
//
// This is the declaration of the VectorWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzPolarizationVector.h>

namespace Herwig {
namespace Helicity {

enum VectorPhase {vector_phase, vector_nophase, default_vector_phase=vector_nophase};

using namespace ThePEG;
using ThePEG::Helicity::LorentzPolarizationVector;

/** \ingroup Helicity
 *  \author Peter Richardson
 *  The VectorWaveFunction class is designed to store the wavefunction
 *  of a vector in a form suitable for use in helicity amplitude calculations 
 *  of the matrix element using a similar philosophy to the FORTRAN HELAS code.
 *
 *  In addition to storing the vector using the LorentzPolarizationVector class
 *  it inherits from the WaveFunctionBase class to provide storage of the 
 *  momentum and particleData for the vector boson.
 *
 *  This class also contains the code which does the actually calculation of the
 *  vector wavefunction.
 *
 *  There are two choices available for the calculation of the wavefunction.
 *  These are set using the VectorPhase enumeration which specifies a default choice.
 *  The first choice, vector_phase, includes a phase factor exp(+/- i phi) for the 
 *  +/- helicity states while the second, vector_nophase, does not.
 *
 *  @see WaveFunctionBase
 *  @see LorentzPolarizationVector
 */
class VectorWaveFunction : public WaveFunctionBase {

public:

  /**
   * Default constructors (set the momentum and Wavefunction)
   */

  /**
   * use a 5-momentum and a LorentzPolarizationVector
   */
  inline VectorWaveFunction(const Lorentz5Momentum &,tcPDPtr,
			    const LorentzPolarizationVector &,Direction=intermediate);

  /**
   * use a 5-momentum and specify all components.
   */
  inline VectorWaveFunction(const Lorentz5Momentum &,tcPDPtr,const Complex &,
			    const Complex &,const Complex &, const Complex &);
  
  /**
   * Use a 5-momentum (specify phase choice).
   */
  inline VectorWaveFunction(const Lorentz5Momentum &,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  /**
   * Set all components of momentum (specify phase choice).
   */
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    int,Direction,VectorPhase=default_vector_phase);

  /**
   * Set 4-momentum components (specify phase choice).
   */
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  /**
   * Set 4-momentum (specify phase choice).
   */
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  /**
   * Set mass zero momentum (specify phase choice).
   */
  inline VectorWaveFunction(Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  /**
   * Set 4 momentum and mass (specify phase choice).
   */
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,int,Direction,
			    VectorPhase=default_vector_phase);

  /**
   * Default constructors (set the momentum and zero the Wavefunction).
   */

  /**
   * Use 5 momentum. 
   */
  inline VectorWaveFunction(Lorentz5Momentum,const tcPDPtr &,Direction); 

  /**
   * Set all components of momentum.
   */
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,Energy,const tcPDPtr &,
			    Direction);

  /**
   * Set 4-momentum components.
   */
  inline VectorWaveFunction(Energy,Energy,Energy,Energy,const tcPDPtr &,Direction);

  /**
   * Set 4-momentum.
   */
  inline VectorWaveFunction(LorentzVector,const tcPDPtr &,Direction);

  /**
   * Set mass zero momentum.
   */
  inline VectorWaveFunction(Energy,const tcPDPtr &,Direction);

  /**
   * Set 4 momentum and mass.
   */
  inline VectorWaveFunction(LorentzVector,Energy,const tcPDPtr &,Direction);

  /**
   * Default constructor.
   */
  inline VectorWaveFunction();

  /**
   * Destructor.
   */
  inline ~VectorWaveFunction();

  /**
   * Subscript operator for the wavefunction.
   */
  inline Complex operator ()(int ) const;

  /**
   * Set components by index.
   */
  inline Complex & operator () (int);

  /**
   * Assignment. 
   */
  inline VectorWaveFunction & operator = (const VectorWaveFunction &);

  /**
   * Return wavefunction as polarization vector. 
   */
  inline LorentzPolarizationVector Wave() const;
  
  /**
   * Get position and time.
   */
  inline Complex x() const;
  inline Complex y() const;
  inline Complex z() const;
  inline Complex t() const;
  
  /**
   * Set position and time.
   */
  inline void setX(const Complex&);
  inline void setY(const Complex&);
  inline void setZ(const Complex&);
  inline void setT(const Complex&);

  /**
   * Reset functions.
   */

  /**
   * Reset the momentum, particle type and direction.
   */
  inline void reset(const Lorentz5Momentum &, const tcPDPtr &, Direction);

  /**
   * Reset the momentum and direction.
   */
  inline void reset(const Lorentz5Momentum &,Direction);

  /**
   * Reset the momentum.
   */
  inline void reset(const Lorentz5Momentum &);

  /**
   * Reset the helicity (recalculation the polarization vector).
   */
  inline void reset(int,VectorPhase=default_vector_phase);

  /**
   * Reset the particle type and direction.
   */
  inline void reset(const tcPDPtr &,Direction);

  /**
   * Reset the particle type.
   */
  inline void reset(const tcPDPtr &);

private:

  /** 
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction();

  /** 
   * Calcuate the wavefunction default phase choice.
   */
  inline void calculateWaveFunction(int);
  
  /**
   * Specify phase choice.
   */
  void calculateWaveFunction(int,VectorPhase=default_vector_phase);

  /**
   * Check particle spin and set pointer.
   */
  inline void checkParticle(const tcPDPtr &);

private:
  
  /**
   * Storage of the wavefunction as a Lorentz Vector.
   */
  LorentzPolarizationVector _wf;
  
};

}
}

#include "VectorWaveFunction.icc"

#endif
