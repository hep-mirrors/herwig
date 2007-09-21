// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPBMultiplicity class.
//

#include "LEPBMultiplicity.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Utilities/StandardSelectors.h"

using namespace Herwig;
using namespace ThePEG;

inline LEPBMultiplicity::LEPBMultiplicity() {
  // B+
  _data[521 ] = MultiplicityInfo(0.175, 0.0044, other);
  // B_s
  _data[531 ] = MultiplicityInfo(0.0450, 0.0038, other);
  // baryons
  _data[5122] = MultiplicityInfo(0.0432, 0.0074, other);
}

void LEPBMultiplicity::analyze(tEventPtr event, long , int , int ) {
  // extract the weakly decaying B hadrons using set to avoid double counting
  set<PPtr> particles;
  map <long,long> eventcount;
  event->select(inserter(particles),WeakBHadronSelector());
  for(set<PPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = abs( (*it)->id() );
    //special for b baryons
    if(ID!=511&&ID!=521&&ID!=531) ID=5122;
    if (_data.find(ID) != _data.end()) {
      eventcount.insert(make_pair(ID,0));
      ++eventcount[ID];
    }
  }
  for(map<long,long>::const_iterator it = eventcount.begin();
      it != eventcount.end(); ++it) {
    _data[it->first].actualCount += it->second;
    _data[it->first].sumofsquares += sqr(double(it->second));
  }
}

NoPIOClassDescription<LEPBMultiplicity> LEPBMultiplicity::initLEPBMultiplicity;
// Definition of the static class description member.

void LEPBMultiplicity::Init() {

  static ClassDocumentation<LEPBMultiplicity> documentation
    ("There is no documentation for the LEPBMultiplicity class");

}

void LEPBMultiplicity::dofinish() {
  string filename = generator()->filename() + ".Bmult";
  ofstream outfile(filename.c_str());
  outfile << 
    "\nParticle multiplicities (compared to LEP data):\n"
    "  ID       Name    simMult     obsMult       obsErr     Sigma\n";
  for (map<long,MultiplicityInfo>::const_iterator it = _data.begin();
       it != _data.end();
       ++it)  {
    MultiplicityInfo multiplicity = it->second;
    string name = (it->first==5122 ? "b baryons" : 
		   generator()->getParticleData(it->first)->PDGName() );
    long N = generator()->currentEventNumber() - 1;
    
    ios::fmtflags oldFlags = outfile.flags();
    outfile << std::scientific << std::showpoint
	    << std::setprecision(3)
	    << setw(7) << it->first << ' '
	    << setw(9) << name << ' ' 
	    << setw(2) << multiplicity.simMultiplicity(N) << " | " 
	    << setw(2) << multiplicity.obsMultiplicity << " +/- " 
	    << setw(2) << multiplicity.obsError << ' '
	    << std::showpos << std::setprecision(1)
	    << multiplicity.nSigma(N) << ' ' 
	    << multiplicity.bargraph(N)
	    << std::noshowpos;
    outfile << '\n';
    outfile.flags(oldFlags);
  }
}
