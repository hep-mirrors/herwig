// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanDalitzAnalysis class.
//

#include "DrellYanDalitzAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig++/Shower2/Kinematics/ShowerParticle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DrellYanDalitzAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DrellYanDalitzAnalysis::~DrellYanDalitzAnalysis() {}

void DrellYanDalitzAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  if(_nout>50000) return;
  tPVector final=event->primaryCollision()->step(1)->getFinalState();
  ParticleVector jets[2];
  for(unsigned int ix=0;ix<final.size();++ix)
    {
      PPtr part=final[ix];
      if(!part->dataPtr()->coloured()) continue;
      do
	{
	  part=part->parents()[0];
	}
      while(dynamic_ptr_cast<ShowerParticlePtr>(part));
      if(part==event->incoming().first) jets[0].push_back(final[ix]);
      else if(part==event->incoming().second) jets[1].push_back(final[ix]);
    }
  Lorentz5Momentum pb,pc,pg;
  ParticleVector::const_iterator cit;
  for(cit=event->incoming().first->children().begin();
      cit!=event->incoming().first->children().end();++cit)
    if((*cit)->dataPtr()->coloured()) pb=(*cit)->momentum();
  for(cit=event->incoming().second->children().begin();
      cit!=event->incoming().second->children().end();++cit)
    if((*cit)->dataPtr()->coloured()) pc=(*cit)->momentum();
  Lorentz5Momentum pjets[2];
  for(unsigned int ix=0;ix<jets[0].size();++ix) pjets[0]+=jets[0][ix]->momentum();
  for(unsigned int ix=0;ix<jets[1].size();++ix) pjets[1]+=jets[1][ix]->momentum();
  // work out the momenta
  if(jets[0].empty()&&jets[1].empty())
    return;
  else if(jets[0].empty())
    {
      pg=pjets[1];
      //      pc-=pg;
    }
  else if(jets[1].empty())
    {
      pg=pjets[0];
      //      pb-=pg;
    }
  else
    {
      if(pjets[0].vect().perp2()>pjets[1].vect().perp2())
	{
	  pg=pjets[0];
	  //	  pb-=pg;
	}
      else
	{
	  pg=pjets[1];
	  //	  pc-=pg;
	}
    }
  Energy2 Q2=(pb+pc-pg).m2();
  Energy2 sbar=(pb+pc).m2()/Q2;
  Energy2 tbar=(pb-pg).m2()/Q2;
  //Energy2 ubar=(pc-pg).m2()/Q2;
  _output[0].push_back(make_pair(sbar,tbar));
}

LorentzRotation DrellYanDalitzAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DrellYanDalitzAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void DrellYanDalitzAnalysis::analyze(tPPtr) {}

void DrellYanDalitzAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DrellYanDalitzAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<DrellYanDalitzAnalysis> DrellYanDalitzAnalysis::initDrellYanDalitzAnalysis;
// Definition of the static class description member.

void DrellYanDalitzAnalysis::Init() {

  static ClassDocumentation<DrellYanDalitzAnalysis> documentation
    ("There is no documentation for the DrellYanDalitzAnalysis class");

}


void DrellYanDalitzAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  ofstream file;
  file.open("DrellYan.top");
  file << "SET WINDOW X 2 9 Y 2 9\n";
  file << "SET LIMITS X 1 10 Y 0 -10\n";
  file << "SET FONT DUPLEX\n";
  file << "TITLE BOTTOM \"sbar\"\n";
  file << "TITLE LEFT \"tbar\"\n";
  for(unsigned int ix=0;ix<_output[0].size();++ix)
    {file << _output[0][ix].first << " " <<  _output[0][ix].second << "\n";}
  file << "PLOT\n";
  //for(unsigned int ix=0;ix<_output[1].size();++ix)
  //  {file << _output[1][ix].first << " " <<  _output[1][ix].second << "\n";}
  //file << "PLOT RED\n";
  // plot the limits
  double kappa[2]={1.,1.};
  double shat,that;
  for(shat=1.;shat<10.;shat+=0.01)
    {
      that=kappa[0]*(1.-shat)/(kappa[0]+shat);
      file << shat << " " << that << "\n";
    }
  file << "join red\n"; 
  for(shat=1.;shat<10.;shat+=0.01)
    {
      that=shat*(1.-shat)/(kappa[1]+shat);
      file << shat << " " << that << "\n";
    }
  file << "join blue\n";
  for(shat=1.;shat<10.;shat+=0.01)
    {
      that=(1.-shat);
      file << shat << " " << that << "\n";
    }
  file << "join\n";
  file << "NEW FRAME\n";
  file << "SET WINDOW X 2 9 Y 2 9\n";
  file << "SET LIMITS X 1 50 Y 0 -50\n";
  file << "SET FONT DUPLEX\n";
  file << "TITLE BOTTOM \"sbar\"\n";
  file << "TITLE LEFT \"tbar\"\n";
  for(unsigned int ix=0;ix<_output[0].size();++ix)
    {file << _output[0][ix].first << " " <<  _output[0][ix].second << "\n";}
  file << "PLOT\n";
  //for(unsigned int ix=0;ix<_output[1].size();++ix)
  //  {file << _output[1][ix].first << " " <<  _output[1][ix].second << "\n";}
  //file << "PLOT RED\n";
  // plot the limits
  for(shat=1.;shat<50.;shat+=0.01)
    {
      that=kappa[0]*(1.-shat)/(kappa[0]+shat);
      file << shat << " " << that << "\n";
    }
  file << "join red\n"; 
  for(shat=1.;shat<50.;shat+=0.01)
    {
      that=shat*(1.-shat)/(kappa[1]+shat);
      file << shat << " " << that << "\n";
    }
  file << "join blue\n";
  for(shat=1.;shat<50.;shat+=0.01)
    {
      that=(1.-shat);
      file << shat << " " << that << "\n";
    }
  file << "join\n";
  file.close();
}
