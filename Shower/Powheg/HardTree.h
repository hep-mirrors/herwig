// -*- C++ -*-
#ifndef HERWIG_HardTree_H
#define HERWIG_HardTree_H
//
// This is the declaration of the HardTree class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/SudakovFormFactor.h"
#include "HardTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The HardTree class is designed to contain the information required
 * to implement the POWHEG approach for Monte Carlo at next-to-leading order.
 */
class HardTree : public Base {

  /**
   *  Output operator for testing
   */
  friend ostream & operator << (ostream &, const HardTree & );

public:

  /**
   * The default constructor.
   */
  HardTree(vector<HardBranchingPtr>,vector<HardBranchingPtr>);

  /**
   *  Match particles in the ShowerTree to branchings in the HardTree
   */
  inline bool connect(ShowerParticlePtr particle, HardBranchingPtr branching) {
    if(_branchings.find(branching)==_branchings.end()) return false;
    _particles[particle]=branching;
    return true;
  }

  /**
   *  Match the prticles in the ShowerTree to the branchings in the HardTree
   */
  bool connect(ShowerTreePtr);
  
  /**
   *  Access the map between the ShowerParticle and the HardBranching
   */
  inline map<ShowerParticlePtr,tHardBranchingPtr> & particles() 
  {return _particles;}

  /**
   *  Access the set of branchings
   */
  inline set<HardBranchingPtr> & branchings() {return _branchings;}
  
  /**
   * Access the incoming branchings
   */
  inline set<HardBranchingPtr> & incoming() {return _spacelike;}

private:

  /**
   *  The ShowerTree
   */
  ShowerTreePtr _tree;

  /**
   *  Map from the particles in the ShowerTree to the HardBranchings
   */
  map<ShowerParticlePtr,tHardBranchingPtr> _particles;

  /**
   *  The HardBranchings in the hard process
   */
  set<HardBranchingPtr> _branchings;

  /**
   *  The HardBranchings which initiate the space-like showers
   */
  set<HardBranchingPtr> _spacelike;

};

/**
 * The HardBranching class is designed to contain the information needed for
 * an individual branching in the POWHEG approach
 */
class HardBranching : public Base {

  /**
   *  The HardTree is friend
   */
  friend class HardTree;

public:

  /**
   * The default constructor
   * @param particle The particle which is branching
   * @param sudakov  The Sudakov form factor for the branching
   * @param parent   The parent for the branching
   * @param incoming Whether the particle is incoming or outgoing
   */
  HardBranching(ShowerParticlePtr particle, SudakovPtr sudakov,
		tHardBranchingPtr parent,bool incoming);

  /**
   * Add a child of the branching
   * @param child The child of the branching
   */
  inline void addChild(HardBranchingPtr child) {_children.push_back(child);}

  /**
   *  Set the momenta of the particles
   */
  void setMomenta(LorentzRotation R, double alpha, Lorentz5Momentum pt,
		  bool setMomentum=true);

  /**
   *  Use the Sudakov to fix the colours
   */
  void fixColours();

  /**
   *  Set and get members for the private member variables
   */
  //@{
  /**
   * Return the branching particle.
   */
  inline tShowerParticlePtr branchingParticle() const {return _particle;}

  /**
   * Set the branching particle
   */
  inline void branchingParticle(ShowerParticlePtr in) {_particle=in;}

  /**
   * Get the original momentum
   */
  inline const Lorentz5Momentum & original() const {return _original;}

  /**
   * Set the original momentum
   */
  inline void original(const Lorentz5Momentum & in) {_original=in;}

  /**
   *  Get the p reference vector
   */
  inline const Lorentz5Momentum & pVector() const {return _p;}

  /**
   *  Set the p reference vector
   */
  inline void pVector(const Lorentz5Momentum & in) {_p=in;}

  /**
   *  Get the n reference vector
   */
  inline const Lorentz5Momentum & nVector() const {return _n;}

  /**
   *  Set the n reference vector
   */
  inline void nVector(const Lorentz5Momentum & in) {_n=in;}

  /**
   *  Get the transverse momentum vector
   */
  inline const Lorentz5Momentum & qPerp() const {return _qt;}

  /**
   *  Set the transverse momentum vector
   */
  inline void qPerp(const Lorentz5Momentum & in) {_qt=in;}

  /**
   *  Get the momentum the particle should have as the start of a shower
   */
  inline const Lorentz5Momentum & showerMomentum() const {return _shower;}

  /**
   *  Get the momentum the particle should have as the start of a shower
   */
  inline void showerMomentum(const Lorentz5Momentum & in ) {_shower=in;}

  /**
   *  Get the transverse momentum
   */
  inline Energy pT() const {return _pt;}

  /**
   *  Set the transverse momentum
   */
  inline void pT(Energy in) { _pt=in;}

  /**
   *  Get whether the branching is incoming or outgoing
   */
  inline bool incoming() const {return _incoming;}

  /**
   *  Set whether the branching is incoming or outgoing
   */
  inline void incoming(bool in) {_incoming=in;}

  /**
   *  The parent of the branching
   */
  inline tHardBranchingPtr parent() const {return _parent;}

  /**
   *  The parent of the branching
   */
  inline void parent(tHardBranchingPtr in) {_parent=in;}

  /**
   *  The Sudakov form-factor
   */
  inline SudakovPtr sudakov() const {return _sudakov;}

  /**
   *  The Sudakov form-factor
   */
  inline void sudakov(SudakovPtr in) {_sudakov=in;}

  /**
   *  Get the beam particle
   */
  inline PPtr beam() const {return _beam;}

  /**
   *  Set the beam particle
   */
  inline void beam(PPtr in) {_beam=in;}

  /**
   * The children
   */
  inline vector<HardBranchingPtr> & children() {return _children;}
  //@}

  /**
   *  Information on the Shower variables for the branching
   */
  //@{
  /**
   *  Get the evolution scale
   */
  inline Energy scale() const {return _scale;}

  /**
   *  The evolution scale
   */
  inline void scale(Energy in) {_scale=in;}

  /**
   *  The energy fraction
   */
  inline double z() const {return _z;}

  /**
   *  The energy fraction
   */
  inline void z(double in) {_z=in;}

  /**
   *  The azimthual angle
   */
  inline double phi() const {return _phi;}

  /**
   *  The azimthual angle
   */
  inline void phi(double in) {_phi=in;}
  //@}

  /**
   *  Colour partners
   */
  //@{
  /**
   *  Get the colour partner
   */
  inline tHardBranchingPtr colourPartner() const {return _partner;}

  /**
   *  The colour partner of the branching
   */
  inline void colourPartner(tHardBranchingPtr in) {_partner=in;}

  //@}

private:

  /**
   *  The branching particle
   */
  ShowerParticlePtr _particle;

  /**
   *  The rescaled momentum
   */
  Lorentz5Momentum _original;

  /**
   *  The \f$p\f$ reference vector
   */
  Lorentz5Momentum _p;

  /**
   *  The \f$n\f$ reference vector
   */
  Lorentz5Momentum _n;

  /**
   *  The transverse momentum vector
   */
  Lorentz5Momentum _qt;

  /**
   *  The momentum the particle should have as the start of a shower
   */
  Lorentz5Momentum _shower;

  /**
   *  The transverse momentum
   */
  Energy _pt;

  /**
   *  Whether the branching is incoming or outgoing
   */
  bool _incoming;

  /**
   *  Information on the Shower variables for the branching
   */
  //@{
  /**
   *  The evolution scale
   */
  Energy _scale;

  /**
   *  The energy fraction
   */
  double _z;

  /**
   *  The azimthual angle
   */
  double _phi;
  //@}

  /**
   *  The parent of the branching
   */
  tHardBranchingPtr _parent;

  /**
   *  The Sudakov form-factor
   */
  SudakovPtr _sudakov;

  /**
   * The children
   */
  vector<HardBranchingPtr> _children;

  /**
   *  The beam particle
   */
  PPtr _beam;

  /**
   *  The colour partner
   */
  tHardBranchingPtr _partner;
};

inline ostream & operator<<(ostream & os, const HardTree & x) {
  os << "Output of HardTree " << &x << "\n";
  for(set<HardBranchingPtr>::const_iterator it=x._branchings.begin();
      it!=x._branchings.end();++it) {
    os << "Hard Particle: " << *(**it).branchingParticle() << " has Sudakov " 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " << (**it).branchingParticle()->colourLine() << "\t" 
       <<  (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children : " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 <<  (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  for(set<HardBranchingPtr>::const_iterator it=x._spacelike.begin();
      it!=x._spacelike.end();++it) {
    os << "SpaceLike: " << *(**it).branchingParticle() << " has Sudakov" 
       << (**it).sudakov() << "\n";
    os << "It's colour lines are " 
       << (**it).branchingParticle()->colourLine() << "\t" 
       << (**it).branchingParticle()->antiColourLine() << "\n";
    for(unsigned int iy=0;iy<(**it).children().size();++iy) {
      os << "\t Children: " << *(**it).children()[iy]->branchingParticle()
	 << "\n";
      os << "It's colour lines are " 
	 << (**it).children()[iy]->branchingParticle()->colourLine() << "\t" 
	 << (**it).children()[iy]->branchingParticle()->antiColourLine() << "\n";
    }
  }
  return os;
}

}

#endif /* HERWIG_HardTree_H */
