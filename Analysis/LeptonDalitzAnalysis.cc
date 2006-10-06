// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonDalitzAnalysis class.
//

#include "LeptonDalitzAnalysis.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LeptonDalitzAnalysis.tcc"
#endif


using namespace Herwig;

LeptonDalitzAnalysis::~LeptonDalitzAnalysis() {}

void LeptonDalitzAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  // Rotate to CMS, extract final state particles and call analyze(particles).
  AnalysisHandler::analyze(event, ieve, loop, state);
  if(_nout>50000) return;
  tPVector final=event->primaryCollision()->step(1)->getFinalState();
  tPVector quark,anti,gluon;
  for(unsigned int ix=0;ix<final.size();++ix)
    {
      tPPtr part=final[ix],last=part->parents()[0];
      do
	{
	  part=last;
	  if(part->parents().empty()) last=tPPtr();
	  else last=part->parents()[0];
	}
      while(!part->parents().empty()&&
	    dynamic_ptr_cast<ShowerParticlePtr>(last));
      if(part->id()==22||part->id()==23&&final[ix]->previous())
	  part=final[ix];
      if(part->id()>0&&part->id()<=6) quark.push_back(final[ix]);
      else if(part->id()<0&&part->id()>=-6) anti.push_back(final[ix]);
      else if(part->id()==21) gluon.push_back(final[ix]);
    }
  // quark jets
  Lorentz5Momentum pquark[2],panti[2];
  double ycut[2]={0.,0.};
  if(gluon.size()!=0)
    {
      Energy eq(0.),eg(0.),ea(0.);
      for(unsigned int ix=0;ix<quark.size();++ix)
	eq+=quark[ix]->momentum().e();
      for(unsigned int ix=0;ix<anti.size();++ix)
	ea+=anti[ix]->momentum().e();
      for(unsigned int ix=0;ix<gluon.size();++ix)
	eg+=gluon[ix]->momentum().e();
      ++_nout;
      eg+=eq+ea;
      _output[1].push_back(make_pair(2.*eq/eg,2.*ea/eg));
      return;
    }
  else if(quark.size()>1)
    {
      _kint.clearMap();
      KtJet::KtEvent ev = KtJet::KtEvent(_kint.convertToKtVectorList(quark), 1, 1, 1);
      ev.findJetsN(2);
      vector<KtJet::KtLorentzVector> ktjets=ev.getJetsPt();
      ycut[0]=ev.getDMerge(1);
      int nquark[2]={0,0},iq;
      for(int ix=0;ix<ev.getNConstituents();++ix)
	{
	  // find jet
	  iq=ktjets[1].contains(*ev.getConstituents()[ix]);
	  if(quark[ix]->id()<0) --nquark[iq];
	  else if(quark[ix]->id()<=6) ++nquark[iq];
	}
      if(nquark[0]>nquark[1])
	{
	  pquark[0]=ktjets[0];
	  pquark[1]=ktjets[1];
	}
      else
	{
	  pquark[1]=ktjets[0];
	  pquark[0]=ktjets[1];
	}
    }
  // antiquark jets
  if(anti.size()>1)
    {
      _kint.clearMap();
      KtJet::KtEvent ev = KtJet::KtEvent(_kint.convertToKtVectorList(anti), 1, 1, 1);
      ev.findJetsN(2);
      vector<KtJet::KtLorentzVector> ktjets=ev.getJetsPt();
      ycut[1]=ev.getDMerge(1);
      int nanti[2]={0,0},iq;
      for(int ix=0;ix<ev.getNConstituents();++ix)
	{
	  // find jet
	  iq=ktjets[1].contains(*ev.getConstituents()[ix]);
	  if(anti[ix]->id()<0) --nanti[iq];
	  else if(anti[ix]->id()<=6) ++nanti[iq];
	}
      if(nanti[0]<nanti[1])
	{
	  panti[0]=ktjets[0];
	  panti[1]=ktjets[1];
	}
      else
	{
	  panti[1]=ktjets[0];
	  panti[0]=ktjets[1];
	}
    }
  double x[2],sum;
  if(quark.size()==1&&anti.size()==1)
    {
      x[0]=1.0;
      x[1]=1.0;
    }
  else if(quark.size()==1&&anti.size()>1)
    {
      x[0]=quark[0]->momentum().e();
      x[1]=panti[0].e();
      sum=x[0]+x[1]+panti[1].e();
      x[0]*=2./sum;
      x[1]*=2./sum;
    }
  else if(quark.size()>1&&anti.size()==1)
    {
      x[0]=pquark[0].e();
      x[1]=anti[0]->momentum().e();
      sum=x[0]+x[1]+pquark[1].e();
      x[0]*=2./sum;
      x[1]*=2./sum;
    }
  else if(quark.size()>1&&anti.size()>1)
    {
      if(ycut[0]>ycut[1])
	{
	  x[0]=pquark[0].e();
	  x[1]=panti[0].e()+panti[1].e();
	  sum=x[0]+x[1]+pquark[1].e();
	}
      else
	{
	  x[0]=pquark[0].e()+pquark[1].e();
	  x[1]=panti[0].e();
	  sum=x[0]+x[1]+panti[1].e();
	}
      x[0]*=2./sum;
      x[1]*=2./sum;
    }
  else
    {
      x[0]=1.;
      x[1]=1.;
      cerr << "testing fails ?? " << quark.size() << " " << anti.size() << endl;
    }
  ++_nout;
  _output[0].push_back(make_pair(x[0],x[1]));

}

LorentzRotation LeptonDalitzAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LeptonDalitzAnalysis::analyze(const tPVector & particles) {
}

void LeptonDalitzAnalysis::analyze(tPPtr) {}

NoPIOClassDescription<LeptonDalitzAnalysis> LeptonDalitzAnalysis::initLeptonDalitzAnalysis;
// Definition of the static class description member.

void LeptonDalitzAnalysis::Init() {

  static ClassDocumentation<LeptonDalitzAnalysis> documentation
    ("There is no documentation for the LeptonDalitzAnalysis class");

}

void LeptonDalitzAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  ofstream file;
  file.open("e+e-Dalitz.top");
  file << "SET WINDOW X 2 9 Y 2 9\n";
  file << "SET FONT DUPLEX\n";
  file << "SET LIMITS X 0 1 Y 0 1\n";
  file << "TITLE BOTTOM \"X011\"\n";
  file << "CASE         \" X X\"\n";
  file << "TITLE LEFT \"X021\"\n";
  file << "CASE       \" X X\"\n";
  for(unsigned int ix=0;ix<_output[0].size();++ix)
    {file << _output[0][ix].first << " " <<  _output[0][ix].second << "\n";}
  file << "PLOT\n";
  for(unsigned int ix=0;ix<_output[1].size();++ix)
    {file << _output[1][ix].first << " " <<  _output[1][ix].second << "\n";}
  file << "PLOT RED\n";
  // plot the limits
  double kb=1.,kc=1.;
  double xc,xb;
  for(double z=0.0;z<=1.0;z+=0.001)
    {
      xc=1.-z*(1.-z)*kb;
      xb=(2.-xc)*0.5+(z-0.5)*xc;
      file << xb << " " << xc << "\n";
    }
  file << "join red" << "\n";
  for(double z=0.0;z<=1.0;z+=0.001)
    {
      xc=1.-z*(1.-z)*kb;
      xb=(2.-xc)*0.5+(z-0.5)*xc;
      file << xc << " " << xb << "\n";
    }
  file << "join red\n";
  file << 0. << " " << 1. << "\n" << 1. << " " << 0. << "\n" << "join\n";
  file.close();
}
