// -*- C++ -*-
#ifndef Herwig_PhaseSpaceChannel_H
#define Herwig_PhaseSpaceChannel_H
//
// This is the declaration of the PhaseSpaceChannel class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "PhaseSpaceMode.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the PhaseSpaceChannel class.
 */
class PhaseSpaceChannel  {

public:

   friend PersistentOStream;
   friend PersistentIStream;
  
/**
 *  Struct for the intermediates in the phase-space channel
 */
struct PhaseSpaceResonance {

  /**
   *   Enum for the jacobians
   */
  enum Jacobian {BreitWigner};
  
   /**
    *   The constructor
    */
PhaseSpaceResonance() {};
   /**
    *   The constructor
    */
PhaseSpaceResonance(cPDPtr part) :particle(part), mass2(sqr(part->mass())), mWidth(part->mass()*part->width()),
				  jacobian(BreitWigner), children(make_pair(0,0))
{};
  /**
   *  The particle
   */
  cPDPtr particle;

  /**
   *  Mass squared
   */
  Energy2 mass2;

  /**
   *  Mass times width
   */
  Energy2 mWidth;

  /**
   *   Type of jacobian
   */
  Jacobian jacobian;

  /**
   *  The children
   */
  pair<int,int> children;

  /**
   *   The final descendents
   */
  vector<int> descendents;
  
};

public:

/**
 *  Default constructior
 */
PhaseSpaceChannel() {};

/** 
 *  Constructor with incoming particles
 */
PhaseSpaceChannel(tPhaseSpaceModePtr inm) : mode_(inm)
{};   

/**
 * If less than zero indicate that this channel is competed. Otherwise
 * signal the parent of the next added parton.
 */
PhaseSpaceChannel & operator , (int o) {
  if(iAdd_<0) iAdd_ = o;
  else if(o>=0) {
     if(intermediates_[iAdd_].children.first==0)
       intermediates_[iAdd_].children.first  = o;
     else
       intermediates_[iAdd_].children.second = o;
     iAdd_=-1;
  }
  else if(o<0) {
     assert(false);
  }
  return *this;
}

public:

  /**
   *  Initialize the channel
   */
  void init(tPhaseSpaceModePtr mode);

  /**
   *  The weight
   */
  const double & weight() const {return weight_;}
  
  /**
   *  Set the weight
   */
  void  weight(double in) {weight_=in;}
  
public:

  /** 
   * Output operator to allow the structure to be persistently written
   * @param os The output stream
   * @param x The intermediate
   */
  inline friend PersistentOStream & operator<<(PersistentOStream & os, 
					       const PhaseSpaceChannel  & x) {
    os << x.intermediates_;
    return os;
  }

  /** 
   * Input operator to allow persistently written data to be read in
   * @param is The input stream
   * @param x The NBVertex 
   */
  inline friend PersistentIStream & operator>>(PersistentIStream & is,
					       PhaseSpaceChannel & x) {
    is >> x.intermediates_;
    return is;
  }
  
  /**
   *  A friend output operator to allow the channel to be outputted for
   * debugging purposes
   */
  friend ostream & operator<<(ostream & os, const PhaseSpaceChannel & channel) {
    // output of the external particles
    // if(!channel.incoming_.second)
    //   os << "Channel for the decay of " << channel.incoming_.first->PDGName()
    // 	<< " -> ";
    // else 
    //   os << "Channel for " << channel.incoming_.first->PDGName()
    // 	<< ", " << channel.incoming_.second->PDGName()
    // 	<< " -> ";
    //   for(unsigned int ix=1;ix<channel._mode->numberofParticles();++ix)
    //     os << channel._mode->externalParticles(ix)->PDGName() << " ";
    os << endl;
    os << "Proceeds in following steps ";
    for(unsigned int ix=0;ix<channel.intermediates_.size();++ix) {
      os << channel.intermediates_[ix].particle->PDGName() << " -> ";
      //     if(channel._intdau1[ix]>0) {
      //       os << channel._mode->externalParticles(channel._intdau1[ix])->PDGName()  
      // 	 << "(" << channel._intdau1[ix]<< ") ";
      //     }
      //     else {
      //       os << channel._intpart[-channel._intdau1[ix]]->PDGName() 
      // 	 << "(" << channel._intdau1[ix]<< ") ";
      //     }
      //     if(channel._intdau2[ix]>0) {
      //       os << channel._mode->externalParticles(channel._intdau2[ix])->PDGName()  
      // 	 << "(" <<channel._intdau2[ix] << ") ";
      //     }
      //     else{
      //       os << channel._intpart[-channel._intdau2[ix]]->PDGName() 
      // 	 << "(" <<channel._intdau2[ix] << ") ";
      //     }
      os << endl;
    }
    return os;
  }

private:

  /**
   *  Find the external particles which are the children of a given resonance
   */
  void findChildren(const PhaseSpaceResonance & res,
		    vector<int> children) {
    if(res.children.first>0)
      children.push_back(res.children.first);
    else
      findChildren(intermediates_[abs(res.children.first)],children);
    if(res.children.second>0)
      children.push_back(res.children.second);
    else
      findChildren(intermediates_[abs(res.children.second)],children);
  }
  
private:

  /**
   *  Pointer to the phase-space mode
   */
  tPhaseSpaceModePtr mode_;
  
  /**
   *  The intermediates
   */
  vector<PhaseSpaceResonance> intermediates_;
  
  /**
   *   Integer to keep track of what we are adding
   */
  int iAdd_ = -1;
  
  /**
   *  The weight
   */
  double weight_;

};

/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The intermediate
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const PhaseSpaceChannel::PhaseSpaceResonance  & x) {
  os << x.particle << ounit(x.mass2,GeV2) << ounit(x.mWidth,GeV2) << x.children << x.descendents;
  return os;
}
  
/** 
 * Input operator to allow persistently written data to be read in
 * @param is The input stream
 * @param x The NBVertex 
 */
inline PersistentIStream & operator>>(PersistentIStream & is,
				      PhaseSpaceChannel::PhaseSpaceResonance & x) {
  is >> x.particle >> iunit(x.mass2,GeV2) >> iunit(x.mWidth,GeV2) >> x.children >> x.descendents;
  return is;
}

}

#endif /* Herwig_PhaseSpaceChannel_H */
