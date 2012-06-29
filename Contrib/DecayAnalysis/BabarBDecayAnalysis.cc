// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BabarBDecayAnalysis class.
//

#include "BabarBDecayAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

IBPtr BabarBDecayAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr BabarBDecayAnalysis::fullclone() const {
  return new_ptr(*this);
}

void BabarBDecayAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
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
  for(set<tcPPtr>::const_iterator it=particles.begin();
      it!=particles.end();++it) {
    long id = (**it).id();
    if(!(abs(id)==ParticleID::B0||abs(id)==ParticleID::Bplus)) continue;
    map <long,long> count;
    vector<tcPPtr> charmed;
    findCharm(*it,charmed);
    for(unsigned int ix=0;ix<charmed.size();++ix) {  
      long idc = id > 0 ? charmed[ix]->id() : -charmed[ix]->id();
      if(abs(id)==ParticleID::B0) {
	if (_b0data.find(idc) != _b0data.end()) ++count[idc];
      }
      else {
	if (_bplusdata.find(idc) != _bplusdata.end()) ++count[idc];
      }
    }
    if(abs(id)==ParticleID::B0) {
      // store info
      for(map<long,MultiplicityInfo>::iterator it = _b0data.begin();
	  it != _b0data.end(); ++it) {
	long currentcount 
	  = count.find(it->first) == count.end() ? 0 : count[it->first];
	it->second.count += currentcount; 
      }
    }
    else {
      // store info
      for(map<long,MultiplicityInfo>::iterator it = _bplusdata.begin();
	  it != _bplusdata.end(); ++it) {
	long currentcount 
	  = count.find(it->first) == count.end() ? 0 : count[it->first];
	it->second.count += currentcount; 
      }
    }
  }
}

NoPIOClassDescription<BabarBDecayAnalysis> BabarBDecayAnalysis::initBabarBDecayAnalysis;
// Definition of the static class description member.

void BabarBDecayAnalysis::Init() {

  static ClassDocumentation<BabarBDecayAnalysis> documentation
    ("There is no documentation for the BabarBDecayAnalysis class");

}

void BabarBDecayAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string filename = generator()->filename() + ".Bdecaymult";
  string dec="B0->";
  ofstream outfile(filename.c_str());
  outfile << 
    "\nParticle multiplicities (compared to Upsilon(4s) data):\n"
    "  ID       Name            simMult     obsMult       obsErr     Sigma\n";
  for (map<long,MultiplicityInfo>::const_iterator it = _b0data.begin();
       it != _bplusdata.end();
       ++it) {
    if(it==_b0data.end()) {
      it=_bplusdata.begin();
      dec="B+->";
    }
    MultiplicityInfo multiplicity = it->second;
    string name = dec+generator()->getParticleData(it->first)->PDGName();
    ios::fmtflags oldFlags = outfile.flags();
    outfile << std::scientific << std::showpoint
	    << std::setprecision(3)
	    << setw(7) << it->first << ' '
	    << setw(17) << name << ' ' 
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

void BabarBDecayAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // B0 multiplicities
  _b0data[ 421 ] = MultiplicityInfo( 0.081,0.015, other);
  _b0data[-421 ] = MultiplicityInfo( 0.474,0.028, other);
  _b0data[ 411 ] = MultiplicityInfo( 0.023,0.013, other);
  _b0data[-411 ] = MultiplicityInfo( 0.369,0.033, other);
  _b0data[ 431 ] = MultiplicityInfo( 0.079,0.014, other);
  _b0data[-431 ] = MultiplicityInfo( 0.103,0.021, other);
  _b0data[ 4122] = MultiplicityInfo( 0.016,0.011, other);
  _b0data[-4122] = MultiplicityInfo( 0.050,0.021, other);
  // B+ multiplicities
  _bplusdata[ 421 ] = MultiplicityInfo( 0.086,0.007, other);
  _bplusdata[-421 ] = MultiplicityInfo( 0.790,0.040, other);
  _bplusdata[ 411 ] = MultiplicityInfo( 0.025,0.005, other);
  _bplusdata[-411 ] = MultiplicityInfo( 0.099,0.012, other);
  _bplusdata[ 431 ] = MultiplicityInfo( 0.079,0.014, other);
  _bplusdata[-431 ] = MultiplicityInfo( 0.011,0.004, other);
  _bplusdata[ 4122] = MultiplicityInfo( 0.021,0.009, other);
  _bplusdata[-4122] = MultiplicityInfo( 0.028,0.011, other);
}

void BabarBDecayAnalysis::findCharm(tcPPtr parent,vector<tcPPtr> & charmed) {
  if(parent->children().empty()) return;
  for(unsigned int ix=0;ix<parent->children().size();++ix) {
    long id = abs(parent->children()[ix]->id());
    if(id==ParticleID::Dplus   || id==ParticleID::D0 ||
       id==ParticleID::D_splus || id==ParticleID::Lambda_cplus) {
      charmed.push_back(parent->children()[ix]);
    }
    findCharm(parent->children()[ix],charmed);
  }
}
