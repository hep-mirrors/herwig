// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Upsilon4SMultiplicityCount class.
//

#include "Upsilon4SMultiplicityCount.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;

inline IBPtr Upsilon4SMultiplicityCount::clone() const {
  return new_ptr(*this);
}

inline IBPtr Upsilon4SMultiplicityCount::fullclone() const {
  return new_ptr(*this);
}

void Upsilon4SMultiplicityCount::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  set<tcPPtr> particles;
  StepVector steps = event->primaryCollision()->steps();
  if (steps.size() > 2) {
    for ( StepVector::const_iterator it = steps.begin()+2;
	  it != steps.end(); ++it ) {
      (**it).select(inserter(particles), ThePEG::AllSelector());
    }
  }
  if(particles.empty()) return;
  map <long,long> eventcount;
  // get most multiplicities
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID =  abs((*it)->id());
    // sort out neutral kaons
    if(ID==ParticleID::K0) continue;
    if(ID==ParticleID::K_L0||ID==ParticleID::K_S0) ID=ParticleID::K0;
    // otherwise add to count
    if (_data.find(ID) != _data.end()) ++eventcount[ID];
  }
  // final-state charged particles
  particles.clear();
  event->selectFinalState(inserter(particles));
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    if((*it)->dataPtr()->charged()) ++eventcount[0];
  }
  // store info
  for(map<long,MultiplicityInfo>::iterator it = _data.begin();
      it != _data.end(); ++it) {
    long currentcount 
      = eventcount.find(it->first) == eventcount.end() ? 0
      : eventcount[it->first];
    it->second.count += currentcount; 
  }
}

NoPIOClassDescription<Upsilon4SMultiplicityCount> Upsilon4SMultiplicityCount::initUpsilon4SMultiplicityCount;
// Definition of the static class description member.

void Upsilon4SMultiplicityCount::Init() {

  static ClassDocumentation<Upsilon4SMultiplicityCount> documentation
    ("There is no documentation for the Upsilon4SMultiplicityCount class");

}

void Upsilon4SMultiplicityCount::dofinish() {
  AnalysisHandler::dofinish();
  string filename = generator()->filename() + ".mult";
  ofstream outfile(filename.c_str());outfile << 
    "\nParticle multiplicities (compared to Upsilon(4s) data):\n"
    "  ID       Name    simMult     obsMult       obsErr     Sigma\n";
  for (map<long,MultiplicityInfo>::const_iterator it = _data.begin();
       it != _data.end();
       ++it) {
    MultiplicityInfo multiplicity = it->second;
    string name = it->first==0 ? "All chgd" : 
      generator()->getParticleData(it->first)->PDGName();
    ios::fmtflags oldFlags = outfile.flags();
    outfile << std::scientific << std::showpoint
	    << std::setprecision(3)
	    << setw(7) << it->first << ' '
	    << setw(9) << name << ' ' 
	    << setw(2) << multiplicity.simMultiplicity() << " | " 
	    << setw(2) << multiplicity.obsMultiplicity << " +/- " 
	    << setw(2) << multiplicity.obsError << ' '
	    << std::showpos << std::setprecision(1)
	    << multiplicity.nSigma() << ' ' 
	    << multiplicity.bargraph()
	    << std::noshowpos;
    outfile << '\n';
    outfile.flags(oldFlags);
  }
  outfile.close();
}

void Upsilon4SMultiplicityCount::doinitrun() {
  AnalysisHandler::doinitrun();
  // all charged
  _data[0] = MultiplicityInfo(10.71,0.18, lightMeson);
  // kaons
  _data[321]  = MultiplicityInfo( 2.*0.789  , 2.*0.025  , strangeMeson);
  _data[311]  = MultiplicityInfo( 2.*0.64   , 2.*0.04   , strangeMeson);
  _data[323]  = MultiplicityInfo( 2.*0.18   , 2.*0.06   , strangeMeson);
  _data[313]  = MultiplicityInfo( 2.*0.146  , 2.*0.026  , strangeMeson);
  // light mesons
  _data[211]  = MultiplicityInfo( 2.*3.58   , 2.*0.07   , lightMeson);
  _data[111]  = MultiplicityInfo( 2.*2.35   , 2.*0.11   , lightMeson);
  _data[221]  = MultiplicityInfo( 2.*0.176  , 2.*0.016  , lightMeson);
  _data[113]  = MultiplicityInfo( 2.*0.21   , 2.*0.05   , lightMeson);
  _data[333]  = MultiplicityInfo( 2.*0.0342 , 2.*0.0013 , lightMeson);
  // baryons
  _data[2212] = MultiplicityInfo( 2.*0.08   , 2.*0.004  , lightBaryon);
  _data[3122] = MultiplicityInfo( 2.*0.04   , 2.*0.005  , lightBaryon);
  _data[3312] = MultiplicityInfo( 2.*0.0027 , 2.*0.006  , lightBaryon);
  _data[4122] = MultiplicityInfo( 2.*0.045  , 2.*0.012  , other);
  _data[4112] = MultiplicityInfo( 4.*0.0046 , 4.*0.0024 , other);
  _data[4222] = MultiplicityInfo( 4.*0.0042 , 4.*0.0024 , other);
  // charmonium
  _data[443]    = MultiplicityInfo( 2.*0.01094, 2.*0.00032, other);
  _data[100443] = MultiplicityInfo( 2.*0.00307, 2.*0.00021, other);
  _data[20443]  = MultiplicityInfo( 2.*0.00386, 2.*0.00027, other);
  _data[445]    = MultiplicityInfo( 2.*0.0013 , 2.*0.0004 , other);
  // D+, D0, D_s+
  _data[411] = MultiplicityInfo( 2.*0.228, 2.*0.014, other);
  _data[421] = MultiplicityInfo( 2.*0.637, 2.*0.03 , other);
  _data[431] = MultiplicityInfo( 2.*0.083, 2.*0.008, other);
  _data[413] = MultiplicityInfo( 2.*0.225, 2.*0.015, other);
  _data[423] = MultiplicityInfo( 2.*0.260, 2.*0.027, other);
  _data[433] = MultiplicityInfo( 2.*0.063, 2.*0.006, other);
}
