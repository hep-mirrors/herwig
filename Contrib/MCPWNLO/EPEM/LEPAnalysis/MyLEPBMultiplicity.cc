// -*- C++ -*-
//
// MyLEPBMultiplicity.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MyLEPBMultiplicity class.
//

#include "MyLEPBMultiplicity.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Utilities/StandardSelectors.h"

using namespace Herwig;
using namespace ThePEG;

string BranchingInfo::bargraph(long N, BranchingInfo den)
{
  if (obsBranching == 0.0) return "     ?     ";
  else if (nSigma(N,den) >= 6.0)  return "-----|---->";
  else if (nSigma(N,den) >= 5.0)  return "-----|----*";
  else if (nSigma(N,den) >= 4.0)  return "-----|---*-";
  else if (nSigma(N,den) >= 3.0)  return "-----|--*--";
  else if (nSigma(N,den) >= 2.0)  return "-----|-*---";
  else if (nSigma(N,den) >= 1.0)  return "-----|*----";
  else if (nSigma(N,den) > -1.0)  return "-----*-----";
  else if (nSigma(N,den) > -2.0)  return "----*|-----";
  else if (nSigma(N,den) > -3.0)  return "---*-|-----";
  else if (nSigma(N,den) > -4.0)  return "--*--|-----";
  else if (nSigma(N,den) > -5.0)  return "-*---|-----";
  else if (nSigma(N,den) > -6.0)  return "*----|-----";
  else                            return "<----|-----";
}

inline MyLEPBMultiplicity::MyLEPBMultiplicity() {
  // B+
  _data[521 ] = BranchingInfo(0.403, 0.009);
  // B_s
  _data[531 ] = BranchingInfo(0.103, 0.009);
  // baryons
  _data[5122] = BranchingInfo(0.091, 0.015);
  // b's
  _data[5]    = BranchingInfo(0.   , 0.   );
}

void MyLEPBMultiplicity::analyze(tEventPtr event, long , int , int ) {
  // extract the weakly decaying B hadrons using set to avoid double counting

  string infile = "EPEM.dat";
  ifstream normin;
  string stringin = "";
  double normfac = 1;
  double stringdoub = 0;
  bool next = 0, filerror = 0;
  normin.open(infile.c_str());
  if(!normin) { cerr << "Error: Failed to open file " << infile << endl; filerror = 1; normfac = 1.; }
  if(!filerror) {
    while(normin) { 
      
      normin >> stringin;
      
      if(next) {   normin >> stringdoub; normfac = stringdoub; break; }
      if(stringin == "NORMFACTOR") { next = 1; }
        
    }
  }
  double eventweight = generator()->currentEvent()->weight()*normfac;
  set<PPtr> particles;
  map <long,long> eventcount;
  StepVector steps = event->primaryCollision()->steps();
  steps[0]->select(inserter(particles), ThePEG::AllSelector());
  unsigned int nb=0;
  for(set<PPtr>::const_iterator cit=particles.begin();cit!=particles.end();++cit) {
    if(abs((*cit)->id())==ParticleID::b) ++nb;
  }
  if(nb!=0) eventcount.insert(make_pair(5,nb));


  particles.clear();
  event->select(inserter(particles),WeakBHadronSelector());
  for(set<PPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = abs( (*it)->id() );
    //special for b baryons
    if(ID!=511&&ID!=521&&ID!=531) ID=5122;
    if (_data.find(ID) != _data.end()) {
      eventcount.insert(make_pair(ID,0));
      ++eventcount[ID]*eventweight;
      //  ++eventcount[ID];
    }
  }
  for(map<long,long>::const_iterator it = eventcount.begin();
      it != eventcount.end(); ++it) {
    _data[it->first].actualCount += it->second*eventweight;
    // _data[it->first].actualCount += it->second;
    _data[it->first].sumofsquares += sqr(double(it->second));
  }
}

NoPIOClassDescription<MyLEPBMultiplicity> MyLEPBMultiplicity::initMyLEPBMultiplicity;
// Definition of the static class description member.

void MyLEPBMultiplicity::Init() {

  static ClassDocumentation<MyLEPBMultiplicity> documentation
    ("There is no documentation for the MyLEPBMultiplicity class");

}

void MyLEPBMultiplicity::dofinish() {
  string filename = generator()->filename() + ".Bmult";
  ofstream outfile(filename.c_str());
  outfile << 
    "\nB branching fraction (compared to LEP data):\n"
    "  ID       Name    simMult     obsMult       obsErr     Sigma\n";
  long N = generator()->currentEventNumber() - 1;
  BranchingInfo den = _data[5];
  for (map<long,BranchingInfo>::const_iterator it = _data.begin();
       it != _data.end();
       ++it)  {
    if(it->first==5) continue;
    BranchingInfo multiplicity = it->second;
    string name = (it->first==5122 ? "b baryons" : 
		   generator()->getParticleData(it->first)->PDGName() ) +"     ";
    ios::fmtflags oldFlags = outfile.flags();
    outfile << std::scientific << std::showpoint
	    << std::setprecision(3)
	    << setw(7) << it->first << ' '
	    << setw(9) << name << ' ' 
	    << setw(2) 
	    << multiplicity.simBranching(N,den) << " | " 
	    << setw(2) 
	    << multiplicity.obsBranching 
	    << " +/- " 
	    << setw(2) << multiplicity.obsError << ' '
	    << std::showpos << std::setprecision(1)
	    << multiplicity.nSigma(N,den) << ' ' 
	    << multiplicity.bargraph(N,den)
	    << std::noshowpos;
    outfile << '\n';
    outfile.flags(oldFlags);
  }
}
