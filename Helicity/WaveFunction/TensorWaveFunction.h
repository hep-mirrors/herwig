// -*- C++ -*-
#ifndef HERWIG_TensorWaveFunction_H
#define HERWIG_TensorWaveFunction_H
//
// This is the declaration of the TensorWaveFunction class.
//
#include "WaveFunctionBase.h"
#include "VectorWaveFunction.h"
#include <ThePEG/Helicity/LorentzTensor.h>

namespace Herwig {
namespace Helicity {
using ThePEG::Helicity::LorentzTensor;

/**\ingroup Helicity
 * Definition of the enumerated values of the phase to include in the 
 * calculation of the polarization tensor.
 */
enum TensorPhase 
{
  tensor_phase, /**< Include the phase factor.*/
  tensor_nophase, /**< No phase-factor. */
  default_tensor_phase=tensor_nophase /**< Default option.*/
};

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The TensorWaveFunction class is designed to store the wavefunction
 *  of a tensor in a form suitable for use in helicity amplitude 
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 * 
 *  In addition to storing the tensor using the LorentzTensor class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and particleData for the tensor particle.
 *
 *  This class also contains the code which does the actually 
 *  calculation of the tensor wavefunction.
 *
 *  There are two choices available for the calculation of the 
 *  wavefunction. These are set using the TensorPhase enumeration 
 *  which specifies a default choice.
 *  The first choice, tensor_phase, includes a phase factor 
 *  \f$\exp(\pm i \phi)\f$ for the \f$\pm\f$ helicity states while the second, 
 *  tensor_nophase, does not.
 *
 *  @see WaveFunctionBase
 *  @see LorentzTensor
 *  @see VectorWaveFunction
 */
class TensorWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and the components of the tensor.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param xx The \f$xx\f$ component.
   * @param xy The \f$xy\f$ component.
   * @param xz The \f$xz\f$ component.
   * @param xt The \f$xt\f$ component.
   * @param yx The \f$yx\f$ component.
   * @param yy The \f$yy\f$ component.
   * @param yz The \f$yz\f$ component.
   * @param yt The \f$yt\f$ component.
   * @param zx The \f$zx\f$ component.
   * @param zy The \f$zy\f$ component.
   * @param zz The \f$zz\f$ component.
   * @param zt The \f$zt\f$ component.
   * @param tx The \f$tx\f$ component.
   * @param ty The \f$ty\f$ component.
   * @param tz The \f$tz\f$ component.
   * @param tt The \f$tt\f$ component.
   */
  inline TensorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,
			    Complex xx,Complex xy,Complex xz,Complex xt,Complex yx,
			    Complex yy,Complex yz,Complex yt,Complex zx,Complex zy,
			    Complex zz,Complex zt,Complex tx,Complex ty,Complex tz,
			    Complex tt);

  /**
   * Constructor, set the momentum, helicity, direction and optionally the phase
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(const Lorentz5Momentum & p,const tcPDPtr & part,int ihel,
			    Direction dir,TensorPhase phase=default_tensor_phase);


  /**
   * Constructor, set the momentum components, mass, helicity and direction,
   * optionally the choice of the phase.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m  The mass.
   * @param part The ParticleData pointer.
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part,int ihel,Direction dir,
			    TensorPhase phase=default_tensor_phase);

  /**
   * Constructor, set the momentum components, helicity and direction,
   * optionally the choice of the phase.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(Energy px,Energy py,Energy pz,Energy E,
			    const tcPDPtr & part,int ihel,Direction dir,
			    TensorPhase phase=default_tensor_phase);

  /**
   * Constructor, set the 4-momentum, helicity, direction and optionally the phase
   * @param p The 4-momentum.
   * @param part The ParticleData pointer
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(LorentzMomentum p,const tcPDPtr & part,int ihel,
			    Direction dir, TensorPhase phase=default_tensor_phase);

  /**
   * Constructor, set the mass, zero the momentum and set the helicity, direction 
   * and optionally the phase.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(Energy m,const tcPDPtr & part,int ihel,Direction dir,
			    TensorPhase phase=default_tensor_phase);

  /**
   * Constructor, set the 4-momentum, mass, helicity, direction and optionally the phase
   * @param p The 4-momentum.
   * @param m The mass.
   * @param part The ParticleData pointer
   * @param ihel The helicity.
   * @param dir The direction.
   * @param phase The phase choice.
   */
  inline TensorWaveFunction(LorentzMomentum p,Energy m,const tcPDPtr & part,int ihel,
			    Direction dir,TensorPhase phase=default_tensor_phase);

  /**
   * Constructor, set the 5-momentum and direction, zero the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(Lorentz5Momentum  p,const tcPDPtr & part,Direction dir); 

  /** 
   * Constructor, set the momentum components, mass and direction, zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(Energy px,Energy py,Energy pz,Energy E,Energy m,
			    const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the momentum components and direction, zero the wavefunction.
   * @param px The x component of the momentum.
   * @param py The y component of the momentum.
   * @param pz The z component of the momentum.
   * @param E  The energy.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(Energy px,Energy py,Energy pz,Energy E,const tcPDPtr & part,
			    Direction dir);

  /**
   * Constructor, set the 4-momentum and direction, zero the wavefunction.
   * @param p The 4-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(LorentzMomentum p,const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the mass and direction, zero the momentum and wavefunction.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(Energy m,const tcPDPtr & part,Direction dir);

  /**
   * Constructor, set the 4-momentum, mass and direction, zero the wavefunction.
   * @param p The 4-momentum.
   * @param m The mass.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline TensorWaveFunction(LorentzMomentum p,Energy m,const tcPDPtr & part,
			    Direction dir);

  /** 
   * Default constructor.
   */
  inline TensorWaveFunction();

  /**
   * Destructor.
   */
  inline ~TensorWaveFunction();
  //@}

  /**
   * Assignment. 
   */
  inline TensorWaveFunction & operator = (const TensorWaveFunction &);

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * Subscript operator for the wavefunction.
   */
  inline Complex operator ()(int,int ) const;

  /**
   * Set components by index.
   */
  inline Complex & operator () (int,int);

  /**
   * Return wavefunction as polarization vector.
   */
  inline LorentzTensor Wave() const;

  /**
   * Get the \f$xx\f$ component.
   */
  inline Complex xx() const;

  /**
   * Get the \f$yx\f$ component.
   */
  inline Complex yx() const;

  /**
   * Get the \f$zx\f$ component.
   */
  inline Complex zx() const;

  /**
   * Get the \f$tx\f$ component.
   */
  inline Complex tx() const;

  /**
   * Get the \f$xy\f$ component.
   */
  inline Complex xy() const;

  /**
   * Get the \f$yy\f$ component.
   */
  inline Complex yy() const;

  /**
   * Get the \f$zy\f$ component.
   */
  inline Complex zy() const;

  /**
   * Get the \f$ty\f$ component.
   */
  inline Complex ty() const;

  /**
   * Get the \f$xz\f$ component.
   */
  inline Complex xz() const;

  /**
   * Get the \f$yz\f$ component.
   */
  inline Complex yz() const;

  /**
   * Get the \f$zz\f$ component.
   */
  inline Complex zz() const;

  /**
   * Get the \f$tz\f$ component.
   */
  inline Complex tz() const;

  /**
   * Get the \f$xt\f$ component.
   */
  inline Complex xt() const;

  /**
   * Get the \f$yt\f$ component.
   */
  inline Complex yt() const;

  /**
   * Get the \f$zt\f$ component.
   */
  inline Complex zt() const;

  /**
   * Get the \f$tt\f$ component.
   */
  inline Complex tt() const;

  /**
   * Set the \f$xx\f$ component.
   */
  inline void setXX(Complex);

  /**
   * Set the \f$yx\f$ component.
   */
  inline void setYX(Complex);

  /**
   * Set the \f$zx\f$ component.
   */
  inline void setZX(Complex);

  /**
   * Set the \f$tx\f$ component.
   */
  inline void setTX(Complex);

  /**
   * Set the \f$xy\f$ component.
   */
  inline void setXY(Complex);

  /**
   * Set the \f$yy\f$ component.
   */
  inline void setYY(Complex);

  /**
   * Set the \f$zy\f$ component.
   */
  inline void setZY(Complex);

  /**
   * Set the \f$ty\f$ component.
   */
  inline void setTY(Complex);

  /**
   * Set the \f$xz\f$ component.
   */
  inline void setXZ(Complex);

  /**
   * Set the \f$yz\f$ component.
   */
  inline void setYZ(Complex);

  /**
   * Set the \f$zz\f$ component.
   */
  inline void setZZ(Complex);

  /**
   * Set the \f$tz\f$ component.
   */
  inline void setTZ(Complex);

  /**
   * Set the \f$xt\f$ component.
   */
  inline void setXT(Complex);

  /**
   * Set the \f$yt\f$ component.
   */
  inline void setYT(Complex);

  /**
   * Set the \f$zt\f$ component.
   */
  inline void setZT(Complex);

  /**
   * Set the \f$tt\f$ component.
   */
  inline void setTT(Complex);
  //@}

  /**
   * Reset functions.
   */
  //@{
  /**
   * Reset the momentum, particle type and direction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline void reset(const Lorentz5Momentum & p, const tcPDPtr & part, Direction dir);

  /** 
   * Reset the momentum and direction
   * @param p The momentum.
   * @param dir The direction
   */
  inline void reset(const Lorentz5Momentum & p,Direction dir);

  /**
   * Reset the momentum.
   * @param p The momentum.
   */
  inline void reset(const Lorentz5Momentum & p);

  /**
   * Reset helicity (recalculate the tensor ).
   * @param ihel The new helicity.
   * @param phase The phase choice.
   */
  inline void reset(int ihel,TensorPhase phase=default_tensor_phase);

  /**
   * Reset particle type and direction.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  inline void reset(const tcPDPtr & part,Direction dir);

  /**
   * Reset particle type.
   * @param part The ParticleData pointer.
   */
  inline void reset(const tcPDPtr & part);
  //@}

private:

  /**
   * Zero the wavefunction.
   */
  inline void zeroWaveFunction();

  /**
   * Calculate the wavefunction.
   * @param ihel The helicity.
   * @param phase The phase choice.
   */
  void calculateWaveFunction(int ihel,TensorPhase phase=default_tensor_phase);

  /**
   * Check particle spin and set pointer.
   * @param part The ParticleData pointer.
   */
  inline void checkParticle(const tcPDPtr & part);

private:

  /**
   * Storage of the wavefunction as a Lorentz Tensor.
   */
  LorentzTensor _wf;

};
}
}

#include "TensorWaveFunction.icc"

#endif
