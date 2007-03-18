// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanHardGenerator class.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "DrellYanHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;
using namespace std;

void DrellYanHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _prefactor << _power;
}

void DrellYanHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _prefactor >> _power;
}

ClassDescription<DrellYanHardGenerator> DrellYanHardGenerator::initDrellYanHardGenerator;
// Definition of the static class description member.

void DrellYanHardGenerator::Init() {

  static ClassDocumentation<DrellYanHardGenerator> documentation
    ("There is no documentation for the DrellYanHardGenerator class");

  static Reference<DrellYanHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &DrellYanHardGenerator::_alphaS, false, false, true, false, false);

}

void DrellYanHardGenerator::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin(); cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    _beams.push_back(cit->first->beam());
    _partons.push_back(cit->first->progenitor()->dataPtr());
  }
  PPtr boson;
  if(tree->outgoingLines().size()==1) {
    boson=tree->outgoingLines().begin()->first->copy();
  }
  else {
    boson=tree->outgoingLines().begin()->first->copy()->parents()[0];
  }
  _yb=0.5*log((boson->momentum().e()+boson->momentum().pz())/
	      (boson->momentum().e()-boson->momentum().pz()));
  _mass=boson->mass();
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
    _yb *=-1.;
  }



  getEvent();

  (*_hyb)+=_yb;
  (*_hplow)+=_pt/GeV;
  (*_hphigh)+=_pt/GeV;
  (*_hyj)+=_yj;

  //cerr<< "hard event (yb,yj,p) = ("<<_yb<<","<<_yj<<","<<_pt/GeV<<")"<<endl;
  
  // at the moment yb (born variables) are generated according to the leading order- this will have to be changed to include virtual emissions- (i.e generate with Bbar)  
}

//this part tests to see if we can handle the event (produced by hard scattering) with the nason implementation

bool DrellYanHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be a quark and an antiquark
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit)
    {part[ix]=cit->first->progenitor();++ix;}
  // check quark and antiquark
  if(!(part[0]->id()>0&&part[0]->id()<6&&part[1]->id()<0&&part[1]->id()>-6)&&
     !(part[1]->id()>0&&part[1]->id()<6&&part[0]->id()<0&&part[0]->id()>-6)) 
    return false;
  // one or two outgoing particles
  if(tree->outgoingLines().size()>2) return false;
  // find the outgoing particles
  ix=0;  
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    {part[ix]=cjt->first->progenitor();++ix;}
  // outgoing particles (1 which is W/Z)
  if(tree->outgoingLines().size()==1&&
     !(part[0]->id()!=ParticleID::gamma||part[0]->id()!=ParticleID::Z0||
       abs(part[0]->id())==ParticleID::Wplus))
    return false;
  else if(tree->outgoingLines().size()==2)
    {
      if(part[0]->parents().empty()||part[1]->parents().empty()) return false;
      if(part[0]->parents()[0]!=part[1]->parents()[0]) return false;
      if(!(part[0]->parents()[0]->id()==ParticleID::gamma||
	   part[0]->parents()[0]->id()==ParticleID::Z0||
	   abs(part[0]->parents()[0]->id())==ParticleID::Wplus)) return false;
    }
  // can handle it
  return true;
}

void DrellYanHardGenerator::doinit() throw(InitException) {
  HardestEmissionGenerator::doinit();
  //cerr << "testing doinit called\n";
}

double DrellYanHardGenerator::getResult() {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 m2(sqr(_mass));
  Energy et=sqrt(m2+sqr(_pt));
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  if(x<0.||x>1.||y<0.||y>1.) return 0.;
  // longitudinal born fractions
  double x1 = _mass*exp( _yb)/sqrt(s);          
  double y1 = _mass*exp(-_yb)/sqrt(s);
  Energy2 th = -sqrt(s)*x*_pt*exp(-_yj);
  Energy2 uh = -sqrt(s)*y*_pt*exp( _yj);
  Energy2 sh = m2-th-uh;
  // pdf part of the weight
  double pdf[4];
  pdf[0]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],m2,x1);
  pdf[1]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],m2,y1);
  pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],sh, x);
  pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],sh, y);
  double wgt;
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    wgt=0.;
  }
  else {
    wgt= pdf[2]*pdf[3]/pdf[0]/pdf[1]*m2/sh;
  }
  // alpha_S part of the weight
  wgt *=_alphaS->ratio(sqr(_pt));
  // matrix element
  wgt *= 2.*(sqr(th-m2)+sqr(uh-m2))*_pt/(sh*uh*th);
  double wgta=4.*_alphaS->value(sqr(_pt))/2./pi*wgt;
  _weightb->addWeighted(_pt/GeV,wgta);
  // divide by the overestimate
  wgt /= _prefactor*pow(GeV/_pt,_power)/GeV;
  _weighta->addWeighted(_pt/GeV,wgt);
  double wgtc=4.*_alphaS->overestimateValue()/2./pi*_prefactor*pow(GeV/_pt,_power)/GeV;
  _weightc->addWeighted(_pt/GeV,wgtc);
  if(wgt>1.) {
    _ptplot.push_back(_pt);
    _yjplot.push_back(_yj-_yb);
//     cerr << "testing violates maximum" << wgt << "\n";
//     cerr << "testing " << _pt << _yb << " " << _yj << "\n";
//     cerr << "testing tree " << x1 << " " << y1 << "\n";
//     cerr << "testing tree " << x  << " " << y  << "\n";
//     cerr << "testing pdf " 
// 	 << pdf[0] << " " << pdf[1] << " " 
// 	 << pdf[2] << " " << pdf[3] << "\n";
  }
  return wgt;
 } 

int DrellYanHardGenerator::getEvent(){
  // pt cut-off
  Energy minp = 0.1*GeV;  
  // maximum pt (half of centre-of-mass energy
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // limits on the rapidity of the boson
  double minyb = -6.0,maxyb = 6.0;
  // multiply the prefactor by other pieces
  double a = 4.*_alphaS->overestimateValue()/6./pi*(maxyj-minyj)*_prefactor/(_power-1.);
  bool reject;
  _pt = maxp;
  // apply the veto algorithm
  do {
    maxp=_pt;
    _pt = GeV/pow(pow(GeV/maxp,_power-1)-log(UseRandom::rnd())/a,1./(_power-1.));
    _yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    double wgt=getResult();
    // compate weight with overestimate
    reject = UseRandom::rnd()>wgt;
    //no emission event if p goes past p min - basically set to outside
    //of the histogram bounds (hopefully hist object just ignores it)
    if(_pt<minp){	
      cerr << "testing veto " << _pt << " " << minp << "\n";
      reject = false;
      _pt=-10;
      _yb=-10;  
      _yj=-10;
    }
  }
  while(reject);
  return 0;
}

void DrellYanHardGenerator::dofinish() {
  //this is called at the end of run - wite histograms to hist.top
  HardestEmissionGenerator::dofinish();

  ofstream hist_out("hist3.top");

  _hyb->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yb"), string(), string(), string(), string("yb"), string());
  _hyj->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yj"), string(), string(), string(), string("yj"), string());
  _hplow->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());
  _hphigh->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());


  _weighta->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("rejection wgt pT"), string(), string(), string(), string("pT/GeV"), string());
  _weightb->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("raw wgt pT"), string(), string(), string(), string("pT/GeV"), string());
  _weightc->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("over wgt pT"), string(), string(), string(), string("pT/GeV"), string());

  hist_out << "NEW FRAME\n";
  for(unsigned int ix=0,N=min(5000,int(_ptplot.size()));ix<N;++ix) {
    hist_out << _yjplot[ix] << " " << _ptplot[ix]/GeV << "\n";
  }
  hist_out << "plot\n";
}

void DrellYanHardGenerator::doinitrun() {
  _hyj= new_ptr(Histogram(-8.0,8.0,100));
  _hplow= new_ptr(Histogram(0.0,5.0,100));
  _hphigh= new_ptr(Histogram(0.0,100.0,100));
  _hyb= new_ptr(Histogram(-6.0,6.0,100));
  _weighta= new_ptr(Histogram(0.0,100.0,100));
  _weightb= new_ptr(Histogram(0.0,100.0,100));
  _weightc= new_ptr(Histogram(0.0,100.0,100));
  //set up histograms in here- this is called at start of run

  HardestEmissionGenerator::doinitrun();
}



