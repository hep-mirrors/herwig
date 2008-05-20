// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LesHouchesOutput class.
//

#include "LesHouchesOutput.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/RemnantParticle.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/RemnantDecayer.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

void LesHouchesParticle::output(ofstream & outfile) {
  ios::fmtflags oldFlags = outfile.flags();
  outfile << std::scientific << std::showpoint
	  << std::setprecision(8)
	  << setw(7) << idup << setw(7) << istup 
	  << setw(7) << mothup[0] << setw(7) << mothup[1] 
	  << setw(7) << icolup[0] << setw(7) << icolup[1] 
	  << setw(16) << pup[0] << setw(16) << pup[1] 
	  << setw(16) << pup[2] << setw(16) << pup[3] 
	  << setw(16) << pup[4] 
	  << std::setprecision(1)
	  << setw(8) << 0. << setw(8) << 0. << "\n";
  outfile.flags(oldFlags);
}

LesHouchesParticle::LesHouchesParticle(tcPPtr in,map<tcColinePtr,unsigned int> &colourMap,
				       unsigned int & nextcolour) 
  : idup(in->id()) {
  pup[0] = in->momentum().x()/GeV;
  pup[1] = in->momentum().y()/GeV;
  pup[2] = in->momentum().z()/GeV;
  pup[3] = in->momentum().e()/GeV;
  pup[4] = in->momentum().mass()/GeV;
  map<tcColinePtr,unsigned int>::const_iterator it;
  if(in->colourLine()) {
    it =colourMap.find(in->colourLine()); 
    if(it!=colourMap.end()) {
      icolup[0] = it->second;
    }
    else {
      icolup[0]=nextcolour;
      colourMap.insert(make_pair(in->colourLine(),nextcolour));
      ++nextcolour;
    }
  }
  else {
    icolup[0]=0;
  }
  if(in->antiColourLine()) {
    it =colourMap.find(in->antiColourLine()); 
    if(it!=colourMap.end()) {
      icolup[1] = it->second;
    }
    else {
      icolup[1]=nextcolour;
      colourMap.insert(make_pair(in->antiColourLine(),nextcolour));
      ++nextcolour;
    }
  }
  else {
    icolup[1]=0;
  }
  mothup[0]=0;
  mothup[1]=0;
}

LesHouchesOutput::LesHouchesOutput(const LesHouchesOutput & x) 
  : AnalysisHandler(x), _filename(x._filename), _mode(x._mode)
{}

void LesHouchesOutput::persistentOutput(PersistentOStream & os) const {
  os << _filename << _mode;
}

void LesHouchesOutput::persistentInput(PersistentIStream & is, int) {
  is >> _filename >> _mode;
}

ClassDescription<LesHouchesOutput> LesHouchesOutput::initLesHouchesOutput;
// Definition of the static class description member.

void LesHouchesOutput::Init() {

  static ClassDocumentation<LesHouchesOutput> documentation
    ("The LesHouchesOutput class is designed to write simple"
     " Les Houches output files");

  static Parameter<LesHouchesOutput,string> interfaceFilename
    ("Filename", "Name of the output file",
     &LesHouchesOutput::_filename, "");

  static Switch<LesHouchesOutput,unsigned int> interfaceMode
    ("Mode",
     "The level of output",
     &LesHouchesOutput::_mode, 0, false, false);
  static SwitchOption interfaceModeNoIntermediates
    (interfaceMode,
     "NoIntermediates",
     "Only include the incoming and outgoing particles",
     0);
  static SwitchOption interfaceModeIntermediates
    (interfaceMode,
     "Intermediates",
     "Add the intermediates as well",
     1);

}

void LesHouchesOutput::dofinish() {
  AnalysisHandler::dofinish();
  // end of the file stuff
  _leshouchesdump << "</LesHouchesEvents>\n";
  _leshouchesdump.close();
}

void LesHouchesOutput::doinitrun() {
  AnalysisHandler::doinitrun();
  // set default filename unless user-specified name exists
  if ( _filename.empty() )
    _filename = generator()->filename() + ".lhe";
  // open the output stream
  _leshouchesdump.open(_filename.c_str()); 
  _leshouchesdump << "<LesHouchesEvents version=\"1.0\">\n";
  _leshouchesdump << "<init>\n";
  EHPtr eh=generator()->eventHandler();
  // beam particles and energies
  _leshouchesdump << eh->incoming().first->id() << "\t"
		  << eh->incoming().second->id() << "\t"
		  << eh->lumiFnPtr()->beamEMaxA()/GeV << "\t"
		  << eh->lumiFnPtr()->beamEMaxB()/GeV << "\t"
    // ignore PDFs
		  << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t"
    // events have unit weight and only 1 process 
		  << 1 << "\t" << 1 << "\n";
  // only one process with arbitary cross section
  _leshouchesdump << 1. << "\t" << "0." << "\t" << "1." << "\t" << 1 << "\n";
  // end of init
  _leshouchesdump << "</init>\n";
}

  
void LesHouchesOutput::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // output the event
  // get the non-remnant incoming particle
  vector<LesHouchesParticle> particles;
  map<tcColinePtr,unsigned int> colourMap;
  map<tcPPtr,unsigned int> particleMap;
  // create the incoming particles
  tPPtr in=event->incoming().first;
  unsigned int nextColour=500;
  do {
    if(in->children().size()>1) {
      for(unsigned int ix=0;ix<in->children().size();++ix) {
	if(!dynamic_ptr_cast<RemPPtr>(in->children()[ix])&&
	   (in->children()[ix]->id()!=ParticleID::gamma&&
	    abs(in->id())==ParticleID::eminus)) {
	  particleMap.insert(make_pair(in->children()[ix],particles.size()));
	  particles.push_back(LesHouchesParticle(in->children()[ix],colourMap,nextColour));
	  particles.back().istup=-1;
	}
      }
    }
    else {
      particleMap.insert(make_pair(in,particles.size()));
      particles.push_back(LesHouchesParticle(in,colourMap,nextColour));
      particles.back().istup=-1;
    }
    if(in==event->incoming().first) in = event->incoming().second;
    else in = tPPtr();
  }
  while(in);
  // and the outgoing particles
  set<tPPtr> outgoing;
  event->selectFinalState(inserter(outgoing));
  for(set<tPPtr>::const_iterator it=outgoing.begin();it!=outgoing.end();++it) {
    if((**it).parents()[0]==event->incoming().first||
       (**it).parents()[0]==event->incoming().second||
       dynamic_ptr_cast<RemPPtr>((**it).parents()[0])) continue;
    particleMap.insert(make_pair(*it,particles.size()));
    particles.push_back(LesHouchesParticle(*it,colourMap,nextColour));
    particles.back().istup = 1;
    if(_mode==0) {
      particles.back().mothup[0] = 1;
      particles.back().mothup[1] = 2;
    }
  }
  if(_mode!=0) {
    for(set<tPPtr>::const_iterator it=outgoing.begin();it!=outgoing.end();++it) {
      map<tcPPtr,unsigned int>::const_iterator jt=particleMap.find(*it);
      if(jt!=particleMap.end()) {
	findIntermediates(*it,colourMap,particleMap,particles,jt->second,nextColour);
      }
    }
    for(unsigned int ix=0;ix<particles.size();++ix)
      if(particles[ix].mothup[0]==0&&particles[ix].istup>0) {
	particles[ix].mothup[0]=1;
	particles[ix].mothup[1]=2;
      }

  }
  // output the event
  Energy2 scalup = event->primarySubProcess()->shat();
  _leshouchesdump << "<event>\n";
  _leshouchesdump << particles.size() << "\t" << 1 << "\t" << 1. << "\t"
		  << sqrt(scalup/GeV2) << "\t" << SM().alphaEM(scalup) << "\t"
		  << SM().alphaS(scalup) << "\n";
  for(unsigned int ix=0;ix<particles.size();++ix) 
    particles[ix].output(_leshouchesdump);
  _leshouchesdump << "</event>\n";
  _leshouchesdump << flush;
}

void LesHouchesOutput::findIntermediates(tPPtr in, map<tcColinePtr,unsigned int> & colourMap,
					 map<tcPPtr,unsigned int> & particleMap,
					 vector<LesHouchesParticle> & particles,
					 unsigned int child,unsigned int & nextColour) {
  map<tcPPtr,unsigned int>::const_iterator it=particleMap.find(in);
  if(it!=particleMap.end()&&particles[it->second].istup<0) return;
  if(it!=particleMap.end()&&particles[it->second].mothup[0]==0&&
     !in->parents().empty()) {
    if(child!=it->second) particles[child].mothup[0]=it->second+1;
    if(particles[it->second].istup>0)
      findIntermediates(in->parents()[0],colourMap,particleMap,
			particles,it->second,nextColour);
    return;
  }
  if(it!=particleMap.end()) {
    if(child!=it->second) particles[child].mothup[0]=it->second+1;
    return;
  }
  else if(in->children().size()>1) {
    if(it!=particleMap.end()) {
      if(particles[it->second].istup>0) 
	particles[child].mothup[0] = it->second+1;
      return;
    }
    else {
      particleMap.insert(make_pair(in,particles.size()));
      particles.push_back(LesHouchesParticle(in,colourMap,nextColour));
      particles.back().istup = 2;
      particles[child].mothup[0] = particles.size();
      findIntermediates(in->parents()[0],colourMap,particleMap,
			particles,particles.size()-1,nextColour);
    }
  }
  else {
    findIntermediates(in->parents()[0],colourMap,particleMap,
		      particles,child,nextColour);
  }
}
