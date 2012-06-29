// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanDalitzAnalysis class.
//

#include "DrellYanDalitzAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DrellYanDalitzAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  if(_nout>50000) return;
  tPVector final=event->primaryCollision()->step(1)->getFinalState();
  ParticleVector jets[3];
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
      else
	{
	  part=part->parents()[0]->children()[0];
	  PPair ptemp=event->primaryCollision()->primarySubProcess()->incoming();


	  if(part==ptemp.first||part==ptemp.second)
	    jets[2].push_back(final[ix]);
	}
    }
  Lorentz5Momentum pb,pc,pg;
  ParticleVector::const_iterator cit;
  for(cit=event->incoming().first->children().begin();
      cit!=event->incoming().first->children().end();++cit)
    if((*cit)->dataPtr()->coloured()) pb=(*cit)->momentum();
  for(cit=event->incoming().second->children().begin();
      cit!=event->incoming().second->children().end();++cit)
    if((*cit)->dataPtr()->coloured()) pc=(*cit)->momentum();
  Lorentz5Momentum pjets[3];
  for(unsigned int iy=0;iy<3;++iy)
    for(unsigned int ix=0;ix<jets[iy].size();++ix) pjets[iy]+=jets[iy][ix]->momentum();
  // work out the momenta
  bool type(false);
  if(!jets[2].empty())
    {
      pg=pjets[2];type=true;
      pb-=pjets[0];
      pc-=pjets[1];
    }
  else if(jets[0].empty()&&jets[1].empty()) return;
  else if(jets[0].empty()) pg=pjets[1];
  else if(jets[1].empty()) pg=pjets[0];
  else
    {
      if(pjets[0].vect().perp2()>pjets[1].vect().perp2())
	{
	  pg=pjets[0];
	  pc-=pjets[1];
	}
      else
	{
	  pg=pjets[1];
	  pb-=pjets[0];
	}
    }
  Energy2 Q2=(pb+pc-pg).m2();
  double sbar=(pb+pc).m2()/Q2;
  double tbar=(pb-pg).m2()/Q2;
  //double ubar=(pc-pg).m2()/Q2;
  ++_nout;
  if(type)
    {
      _output[1].push_back(make_pair(sbar,tbar));
    }
  else
    {
      _output[0].push_back(make_pair(sbar,tbar));
    }
}

LorentzRotation DrellYanDalitzAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DrellYanDalitzAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void DrellYanDalitzAnalysis::analyze(tPPtr) {}

void DrellYanDalitzAnalysis::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DrellYanDalitzAnalysis::persistentInput(PersistentIStream &, int) {
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
  string fname = generator()->filename() + string("-") + name() + string(".top");
  file.open(fname.c_str());
  file << "SET WINDOW X 2 9 Y 2 9\n";
  file << "SET LIMITS X 1 10 Y 0 -10\n";
  file << "SET FONT DUPLEX\n";
  file << "TITLE BOTTOM \"sbar\"\n";
  file << "TITLE LEFT \"tbar\"\n";
  for(unsigned int ix=0;ix<_output[0].size();++ix)
    {file << _output[0][ix].first << " " <<  _output[0][ix].second << "\n";}
  file << "PLOT\n";
  for(unsigned int ix=0;ix<_output[1].size();++ix)
    {file << _output[1][ix].first << " " <<  _output[1][ix].second << "\n";}
  if(!_output[1].empty()) file << "PLOT RED\n";
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
  for(unsigned int ix=0;ix<_output[1].size();++ix)
    {file << _output[1][ix].first << " " <<  _output[1][ix].second << "\n";}
  if(!_output[1].empty()) file << "PLOT RED\n";
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
