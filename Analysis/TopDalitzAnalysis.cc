// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TopDalitzAnalysis class.
//

#include "TopDalitzAnalysis.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDalitzAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower2/Kinematics/ShowerParticle.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

TopDalitzAnalysis::~TopDalitzAnalysis() {}

void TopDalitzAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  if(_nout>50000) return;
  ParticleSet pert=event->primaryCollision()->step(1)->all();
  tPVector final=event->primaryCollision()->step(1)->getFinalState();
  ParticleSet::const_iterator pit;
  for(pit=pert.begin();pit!=pert.end();++pit)
    {
      // must be top
      if(abs((*pit)->id())!=6) continue;
      // must have two children
      if((*pit)->children().size()!=2) continue;
      // neither should be top
      if(abs((*pit)->children()[0]->id())==6||
	 abs((*pit)->children()[1]->id())==6) continue;
      topShower(*pit,final);
    }
}

LorentzRotation TopDalitzAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void TopDalitzAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void TopDalitzAnalysis::analyze(tPPtr) {}

void TopDalitzAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void TopDalitzAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<TopDalitzAnalysis> TopDalitzAnalysis::initTopDalitzAnalysis;
// Definition of the static class description member.

void TopDalitzAnalysis::Init() {

  static ClassDocumentation<TopDalitzAnalysis> documentation
    ("There is no documentation for the TopDalitzAnalysis class");

}


void TopDalitzAnalysis::topShower(PPtr top,tPVector final)
{
  PPtr orig=top->parents()[0];
  // find the original version of this top
  do
    {
      orig=orig->parents()[0];
    }
  while(dynamic_ptr_cast<ShowerParticlePtr>(orig)&&orig->children().size()==2);
  if(abs(orig->parents()[0]->id())!=ParticleID::t) orig=PPtr();
  // find the W
  PPtr Wboson,borig;
  for(unsigned int ix=0;ix<top->children().size();++ix)
    {
      if(abs(top->children()[ix]->id())==ParticleID::Wplus)
	{
	  Wboson=top->children()[ix];
	  if(abs(Wboson->children()[0]->id())==ParticleID::Wplus)
	    Wboson=Wboson->children()[0];
	}
      else
	{
	  borig=top->children()[ix];
	  if(abs(borig->children()[0]->id())==ParticleID::b)
	    borig=borig->children()[0];
	}
    }
  // find the particles
  tPVector tprod,bottom;
  for(unsigned int ix=0;ix<final.size();++ix)
    {
      tPPtr part=final[ix],last=part->parents()[0];
      tShowerParticlePtr shower;
      do
	{
	  part=last;
	  if(part->parents().empty())
	    {
	      last=tPPtr();
	      shower=tShowerParticlePtr();
	    }
	  else if(last&&last==orig)
	    {break;}
	  else
	    {
	      last=part->parents()[0];
	      shower=dynamic_ptr_cast<ShowerParticlePtr>(last);
	    }
	}
      while(!part->parents().empty()&&shower);
      if(last)
	{
	  if(abs(last->id())==ParticleID::b&&last->parents()[0]==top)
	    bottom.push_back(final[ix]);
	  else if(last==orig)
	    tprod.push_back(final[ix]);
	}
    }
  if(bottom.empty()) bottom.push_back(borig);
  for(unsigned int ix=0;ix<tprod.size();++ix)
    {bottom.push_back(tprod[ix]);}
  // check the top
  Lorentz5Momentum ptotal;
  for(unsigned int ix=0;ix<bottom.size();++ix)
    ptotal+=bottom[ix]->momentum();
  ptotal+=Wboson->momentum();
  ptotal.rescaleMass();
  if(bottom.size()>1)
    {
      _kint.clearMap();
      KtJet::KtEvent ev = KtJet::KtEvent(_kint.convertToKtVectorList(bottom), 1, 1, 1);
      ev.findJetsN(2);
      vector<KtJet::KtLorentzVector> ktjets=ev.getJetsPt();
      int nquark[2]={0,0},iq;
      for(int ix=0;ix<ev.getNConstituents();++ix)
	{
	  // find jet
	  iq=ktjets[1].contains(*ev.getConstituents()[ix]);
	  if(bottom[ix]->id()<0) --nquark[iq];
	  else if(bottom[ix]->id()<=6) ++nquark[iq];
	}
      Lorentz5Momentum pb,pg;
      // identify the jets
      if(top->id()>0)
	{
	  if(nquark[0]>nquark[1])
	    {
	      pb=ktjets[0];
	      pg=ktjets[1];
	    }
	  else 
	    {
	      pb=ktjets[1];
	      pg=ktjets[0];
	    }
	}
      else
	{
	  if(nquark[0]>nquark[1])
	    {
	      pb=ktjets[1];
	      pg=ktjets[0];
	    }
	  else 
	    {
	      pb=ktjets[0];
	      pg=ktjets[1];
	    }
	}
      // boost to the rest frame
      Hep3Vector boost;
      if(orig) boost=-orig->momentum().boostVector();
      else boost=-top->momentum().boostVector();
      pg.boost(boost);
      pb.boost(boost);
      Energy mt(top->mass());
      double xg(2.*pg.e()/mt),xb(2.*pb.e()/mt);
      _output << xg << " " << 2.-xb-xg << "\n";
      ++_nout;
    }
}

void TopDalitzAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  _output << "PLOT\n";
  Energy mt=getParticleData(ParticleID::t)->mass();
  Energy mb=getParticleData(ParticleID::b)->constituentMass();
  Energy mw=getParticleData(ParticleID::Wplus)->mass();
  Energy2 mb2(mb*mb),mt2(mt*mt),mw2(mw*mw);
  Energy2 m122(sqr(mb+mw)),step;
  step=(sqr(mt)-m122)/200.;
  vector<double> upper,lower,xgg;
  for(;m122<=sqr(mt);m122+=step)
    {
      Energy m12=sqrt(m122);
      Energy E2s=0.5*(m122-mb2+mw2)/m12;
      Energy E3s=0.5*(mt2-m122)/m12;
      Energy2 m23max=2.*E2s*E3s+mw2+2.*E3s*sqrt(sqr(E2s)-mw2);
      Energy2 m23min=2.*E2s*E3s+mw2-2.*E3s*sqrt(sqr(E2s)-mw2);
      xgg.push_back(1.-m122/mt2);
      upper.push_back((m122+m23max-mb2)/mt2);
      lower.push_back((m122+m23min-mb2)/mt2);
    }
  for(unsigned int ix=0;ix<upper.size();++ix)
    {_output << xgg[ix] << " " << upper[ix] << "\n";}
  for(int ix=lower.size()-1;ix>=0;--ix)
    {_output << xgg[ix] << " " << lower[ix] << "\n";}
  _output << "JOIN RED " << "\n";
  // phase space for radiation from bottom
  double a=mw2/mt2,c=mb2/mt2,xa,xc,r,xg;
  double lam=sqrt(sqr(1.+a-c)-4.*a);
  // maximal b choice
  //double kappa=4.*(1.-c-2.*sqrt(a)+a);
  // symmetric choice
  double kappa=0.5*(1-a+c+lam)+c;
  // smooth choice
  //double kappa=sqrt(c)*lam*(1.+c-a+lam)/(1+c-a+lam-2.*sqrt(c));
  cerr << "\nbottom kappa " << kappa << endl;
  double xgmax=1.-sqr(sqrt(a)+sqrt(c));
  for(double z=sqrt(c/kappa);z<=1.;z+=0.001)
    {
      xa=1.+a-c-z*(1.-z)*kappa;
      r =0.5*(1.+c/(1.+a-xa));
      xc=(2.-xa)*r+(z-r)*sqrt(sqr(xa)-4.*a);
      xg=(2.-xa)*(1.-r)-(z-r)*sqrt(sqr(xa)-4.*a);
      if(xg<xgmax) _output << xg << " " << xa << "\n";
    }
  _output << "JOIN BLUE" << endl;
  // phase space for radiation from top
  kappa=1+0.25*sqr(1.-a+c+lam)/(kappa-c);
  cerr << "top    kappa " << kappa << endl;
  double u,w,v;
  double zmin=1.-(1.-a)/(kappa+2.*sqrt(a*(kappa-1.)))+0.00001;
  for(double z=0.;z<=1.;z+=0.0001)
    {
      double kmax=2*a + (-1 + a + c)/(-1 + z) - 
	(2*sqrt(a*(1 + c + a*(-1 + z) - z)*pow(-1 + z,2)*z))/pow(-1 + z,2);
      if(kmax<kappa)
	{
	  u = 1+a-c-(1.-z)*kmax;
	  w = 1.-(1.-z)*(kmax-1.);
	  v = 0.;
	  xa =0.5*((u+v)/w+(u-v)/z);
	  xc = w+z-xa;
	  xg = (1.-z)*kmax;
	}
      else
	{
	  u = 1+a-c-(1.-z)*kappa;
	  w = 1.-(1.-z)*(kappa-1.);
	  v = sqrt(sqr(u)-4.*a*w*z);
	  xa =0.5*((u+v)/w+(u-v)/z);
	  xc = w+z-xa;
	  xg = (1.-z)*kappa;
        }
      if(xg<xgmax) _output << xg << " " << xa << "\n";
    }
  _output << "JOIN GREEN" << endl;
  _output.close();
}



