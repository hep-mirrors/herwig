// -*- C++ -*-
#ifndef Herwig_OneLoopFFAWZVertex_H
#define Herwig_OneLoopFFAWZVertex_H
//
// This is the declaration of the OneLoopFFAWZVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "OneLoopFFAWZVertex.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The OneLoopFFAWZVertex class provides the one-loop EW renormalised
 * vertex for the coupling of photons, W and Z bosons to fermions
 * in the Standard Model.
 *
 * @see \ref OneLoopFFAWZVertexInterfaces "The interfaces"
 * defined for OneLoopFFAWZVertex.
 */
class OneLoopFFAWZVertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  OneLoopFFAWZVertex();
  
public:

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const SpinorBarWaveFunction & sbar2,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) {
    if(mass.real()<ZERO) {
      long id = abs(out->id());
      if(id==ParticleID::gamma)             mass = ZERO;
      else if(out->id()==ParticleID::Z0)    mass = muZ_;
      else if(out->id()==ParticleID::Wplus) mass = muW_;
      else assert(false);
      width = ZERO;
    }
    return FFVVertex::evaluate(q2,iopt,out,sp1,sbar2,mass,width);
  }

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const VectorWaveFunction & vec3) {
    return FFVVertex::evaluate(q2,sp1,sbar2,vec3);
  }

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
					 const SpinorBarWaveFunction & sbar2,
					 const VectorWaveFunction & vec3,
					 complex<Energy> mass=-GeV,
					 complex<Energy> width=-GeV) {
    return FFVVertex::evaluate(q2,iopt,out,sbar2,vec3,mass,width);
  }

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const VectorWaveFunction & vec3,
				      complex<Energy> mass=-GeV,
				      complex<Energy> width=-GeV) {
    return FFVVertex::evaluate(q2,iopt,out,sp1,vec3,mass,width);
  }

  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

  /**
   *  Whether to return the LO result or vertex correction
   */
  void setOrder(unsigned int order) {
    assert(order<=1);
    order_ = order;
  }

  VectorWaveFunction selfEnergyCorrection(tcPDPtr particle,
					  const VectorWaveFunction & old);

  complex<InvEnergy2> neutralCurrentBT(int hel1, int hel2,
				       complex<Energy2> mV2,
				       complex<Energy2> mVp2, 
				       Energy2 mQ2, Energy2 mq2, Energy2 ml2,
				       Energy2 sHat, Energy2 tHat, Energy2 uHat);
  
  complex<InvEnergy2> neutralCurrentBU(int hel1, int hel2,
				       complex<Energy2> mV2,
				       complex<Energy2> mVp2, 
				       Energy2 mQ2, Energy2 mq2, Energy2 ml2,
				       Energy2 sHat, Energy2 tHat, Energy2 uHat) {
    return -neutralCurrentBT(hel1,-hel2,mV2,mVp2,mQ2,mq2,ml2,sHat,uHat,tHat);
  }

  vector<vector<complex<InvEnergy2> > > 
  neutralCurrentFBox(tcPDPtr q1, tcPDPtr q2,
		     tcPDPtr l1, tcPDPtr l2,
		     Energy2 sHat, Energy2 tHat, Energy2 uHat);

  /**
   *  Test of the matrix element for Drell-Yan
   */
  void neutralCurrentME(tcPDPtr q1, tcPDPtr q2,
			tcPDPtr l1, tcPDPtr l2,
			Energy2 sHat, Energy2 tHat, Energy2 uHat);
  
public:
  
    /**
   * Renormalisated self energies
   */
  //@{
  /**
   * \f$\hat{\Sigma}_{T}^{\gamma\gamma}(k^2)\f$
   */
  complex<Energy2> SigmaHat_T_AA(Energy2 k2) {
    return Sigma_T_AA(k2)+k2*dZ_AA_;
  }

  /**
   * \f$\hat{\Sigma}_{T}^{\gamma Z}(k^2)\f$
   */
  complex<Energy2> SigmaHat_T_AZ(Energy2 k2) {
    return Sigma_T_AZ(k2)+0.5*dZ_AZ_*k2+0.5*dZ_ZA_*(k2-muZ2_);
  }

  /**
   * \f$\hat{\Sigma}_{T}^{ZZ}(k^2)\f$
   */
  complex<Energy2> SigmaHat_T_ZZ(Energy2 k2) {
    return Sigma_T_ZZ(k2)-dmuZ2_+dZ_ZZ_*(k2-muZ2_);
  }

  /**
   * \f$\hat{\Sigma}_{T}^W(k^2)\f$
   */
  complex<Energy2> SigmaHat_T_W(Energy2 k2) {
    return Sigma_T_W(k2)-dmuW2_+dZ_W_*(k2-muW2_);
  }
  //@}

  /**
   *  The unrenormalised self energies of the gauge bosons
   */
  //@{
  /**
   * \f$\Sigma_{T}^{\gamma\gamma}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AA(Energy2 k2) {
    return Sigma_T_AA_Boson(k2)+Sigma_T_AA_Fermion(k2);
  }

  /**
   * \f$\Sigma_{T}^{\gamma Z}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AZ(Energy2 k2) {
    return Sigma_T_AZ_Boson(k2)+Sigma_T_AZ_Fermion(k2);
  }

  /**
   * \f$\Sigma_{T}^{ZZ}(k^2)\f$
   */
  complex<Energy2> Sigma_T_ZZ(Energy2 k2) {
    return Sigma_T_ZZ_Boson(k2)+Sigma_T_ZZ_Fermion(k2);
  }

  /**
   * \f$\Sigma_{T}^W(k^2)\f$
   */
  complex<Energy2> Sigma_T_W(Energy2 k2) {
    return Sigma_T_W_Boson(k2)+Sigma_T_W_Fermion(k2);
  }
  //@}

  /**
   *  The unrenormalised bosonic contributions to the 
   *  self energies of the gauge bosons
   */
  //@{
  /**
   * Bosonic contribution to \f$\Sigma_{T}^{\gamma\gamma}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AA_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\Sigma_{T}^{\gamma Z}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AZ_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\Sigma_{T}^{ZZ}(k^2)\f$
   */
  complex<Energy2> Sigma_T_ZZ_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\Sigma_{T}^W(k^2)\f$
   */
  complex<Energy2> Sigma_T_W_Boson(Energy2 k2);
  //@}

  /**
   *  The unrenormalised fermionic contributions to the 
   *  self energies of the gauge bosons
   */
  //@{
  /**
   * Fermionic contribution to \f$\Sigma_{T}^{\gamma\gamma}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AA_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\Sigma_{T}^{\gamma Z}(k^2)\f$
   */
  complex<Energy2> Sigma_T_AZ_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\Sigma_{T}^{ZZ}(k^2)\f$
   */
  complex<Energy2> Sigma_T_ZZ_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\Sigma_{T}^W(k^2)\f$
   */
  complex<Energy2> Sigma_T_W_Fermion(Energy2 k2);
  //@}

  /**
   * Derivatives of the boson self energies
   */
  //@{
  /**
   * \f$\frac{\partial \Sigma_{T}^{\gamma\gamma}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AA(Energy2 k2) {
    return DSigma_T_AA_Boson(k2)+DSigma_T_AA_Fermion(k2);
  }

  /**
   * \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AZ(Energy2 k2) {
    return DSigma_T_AZ_Boson(k2)+DSigma_T_AZ_Fermion(k2);
  }

  /**
   * \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_ZZ(Energy2 k2) {
    return DSigma_T_ZZ_Boson(k2)+DSigma_T_ZZ_Fermion(k2);
  }

  /**
   * \f$\frac{\partial \Sigma_{T}^W(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_W(Energy2 k2) {
    return DSigma_T_W_Boson(k2)+DSigma_T_W_Fermion(k2);
  }
  //@}

  /**
   * Bosonic contributions to the derivatives of the boson self energies
   */
  //@{
  /**
   * Bosonic contribution to 
   * \f$\frac{\partial \Sigma_{T}^{\gamma\gamma}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AA_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AZ_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_ZZ_Boson(Energy2 k2);

  /**
   * Bosonic contribution to \f$\frac{\partial \Sigma_{T}^W(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_W_Boson(Energy2 k2);
  //@}

  /**
   * Fermionic contributions to the derivatives of the boson self energies
   */
  //@{
  /**
   * Fermionic contribution to 
   * \f$\frac{\partial \Sigma_{T}^{\gamma\gamma}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AA_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_AZ_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\frac{\partial \Sigma_{T}^{ZZ}(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_ZZ_Fermion(Energy2 k2);

  /**
   * Fermionic contribution to \f$\frac{\partial \Sigma_{T}^W(k^2)}{\partial k^2}\f$
   */
  Complex DSigma_T_W_Fermion(Energy2 k2);
  //@}

  /**
   *  Renormalised loop correction of the photon and Z fermion-fermion vertices
   */
  pair<Complex,Complex> renormalisedVertex(Energy2 sHat,tcPDPtr aa,tcPDPtr,tcPDPtr cc);

  /**
   *  Counter term for the evolution of the EM coupling from \f$q^2=0\f$ to
   *  \f$q^2=M^2_Z\f$.
   */
  Complex deltaAlphaMZ();


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

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OneLoopFFAWZVertex & operator=(const OneLoopFFAWZVertex &);

private:

  /**
   *  The order to return
   */
  unsigned int order_;

  /**
   * Parameters for the EW corrections
   */
  //@{
  /**
   *  EW scheme
   */
  unsigned int EWscheme_;

  /**
   *  Z mass
   */
  Energy mZ_;

  /**
   *  W mass
   */
  Energy mW_;

  /**
   *  Complex \f$m_Z^2\f$
   */
  complex<Energy2> muZ2_;

  /**
   *  Complex \f$m_Z\f$
   */
  complex<Energy>  muZ_;

  /**
   *  Complex \f$m_W^2\f$
   */
  complex<Energy2> muW2_;

  /**
   *  Complex \f$m_W^2\f$
   */
  complex<Energy>  muW_;

  /**
   *   Higgs mass 
   */
  Energy mH_;

  /**
   *   Higgs mass 
   */
  Energy2 mH2_;

  /**
   *  Complex \f$\cos^2\theta_W\f$
   */
  Complex cw2_;

  /**
   *  Complex \f$\cos\theta_W\f$
   */
  Complex cw_;

  /**
   *  Complex \f$\sin^2\theta_W\f$
   */
  Complex sw2_;

  /**
   *  Complex \f$\sin\theta_W\f$
   */
  Complex sw_;

  /**
   *  Electric charges
   */ 
  vector<Complex> ef_;

  /**
   *  Left charges
   */ 
  vector<Complex> gl_;

  /**
   *  Right charges
   */ 
  vector<Complex> gr_;

  /**
   *  Fermion masses
   */
  vector<Energy> fermionMasses_;

  /**
   * The fixed value of \f$\alpha_{EM}\f$
   */
  double alphaEW_;

  /**
   *  The electric charge
   */
  double e_;

  /**
   *  \f$\mathcal{Z}_{AA}\f$ counterterm
   */
  Complex dZ_AA_;

  /**
   *  \f$\mathcal{Z}_{ZA}\f$ counterterm
   */
  Complex dZ_ZA_;

  /**
   *  \f$\mathcal{Z}_{WW}\f$ counterterm
   */
  Complex dZ_W_;

  /**
   *  \f$\mathcal{Z}_{AZ}\f$ counterterm
   */
  Complex dZ_AZ_;

  /**
   *  \f$\mathcal{Z}_{ZZ}\f$ counterterm
   */
  Complex dZ_ZZ_;

  /**
   *  \f$\delta\mu^2_W\f$ counterterm 
   */
  complex<Energy2>  dmuW2_;

  /**
   *  \f$\delta\mu^2_Z\f$ counterterm 
   */
  complex<Energy2>  dmuZ2_;

  /**
   *  \f$\frac{\delta e}{e}\f$ counterterm 
   */
  Complex de_;

  /**
   *  \f$\frac{\delta\sin\theta_W}{\sin\theta_W}\f$ counterterm 
   */
  Complex dsw_;
  //@}

};

}

#endif /* Herwig_OneLoopFFAWZVertex_H */
