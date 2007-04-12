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
  os << _alphaS << _power;
}

void DrellYanHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _power;
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


NasonTreePtr DrellYanHardGenerator::generateHardest(ShowerTreePtr tree,EvolverPtr) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  _quarkplus=true;
  for(cit=tree->incomingLines().begin(); cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    _beams.push_back(cit->first->beam());
    _partons.push_back(cit->first->progenitor()->dataPtr());
    if(cit->first->progenitor()->id()>0&&cit->first->progenitor()->momentum().z()<0)
      _quarkplus=false;
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
  _yb *= _quarkplus ? 1. : -1.;
  _mass=boson->mass();
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }

  getEvent();
  
  //  _ptplot.push_back(_pt/GeV);
  //_yjplot.push_back(_yj);
  //_xplot.push_back(_x);
  //_yplot.push_back(_y);
  (*_hyb)+=_yb;
  (*_hplow)+=_pt/GeV;
  (*_hphigh)+=_pt/GeV;
  (*_hyj)+=_yj;
}
   
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
}

double DrellYanHardGenerator::getResult(int emis_type, Energy pt, double yj) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 m2(sqr(_mass));
  Energy2 scale = m2+sqr(pt);
  Energy  et=sqrt(scale);
  // longitudinal real correction fractions
  _x  = pt*exp( yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  _y  = pt*exp(-yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  // reject if outside region
  if(_x<0.||_x>1.||_y<0.||_y>1.||_x*_y<m2/s) return 0.;
  // longitudinal born fractions
  _x1 = _mass*exp( _yb)/sqrt(s);          
  _y1 = _mass*exp(-_yb)/sqrt(s);
  // mandelstam variables
  Energy2 th = -sqrt(s)*_x*pt*exp(-yj);
  Energy2 uh = -sqrt(s)*_y*pt*exp( yj);
  Energy2 sh = m2-th-uh;
  double res;
  // pdf part of the cross section
  double pdf[4];
  pdf[0]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],m2,_x1);
  pdf[1]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],m2,_y1);
  //qqbar2Zg
  if(emis_type==0) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,_x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,_y);
    res = 4./3./pi*(sqr(th-m2)+sqr(uh-m2))*pt/(sh*uh*th)*GeV;
    res=0.;
  }
  //qg2Zq
  else if(emis_type==1) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,_x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],getParticleData(ParticleID::g),scale,_y);
    res = -1./2./pi*(sqr(uh-m2)+sqr(sh-m2))*pt/(sh*sh*uh)*GeV;
    res=0.;
  }
  //qbarg2Zqbar
  else {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],getParticleData(ParticleID::g),scale,_x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,_y);
    res =- 1./2./pi*(sqr(th-m2)+sqr(sh-m2))*pt/(sh*sh*th)*GeV;
  }  
  //deals with pdf zero issue at large x
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    res=0.;
  }
  else {
    res*=pdf[2]*pdf[3]/pdf[0]/pdf[1]*m2/sh;
  }
  res*=_alphaS->ratio(scale);
  if(res*sqr(pt/GeV)>_max[emis_type]) _max[emis_type]=res*sqr(pt/GeV);
  return res;
 } 

int DrellYanHardGenerator::getEvent(){

  // pt cut-off
  Energy minp = 40.*GeV;  
  // maximum pt (half of centre-of-mass energy)
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  //Working Variables
  Energy pt;
  Energy winning_pt =0.0*GeV;
  double yj;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  bool reject;
  double wgt;
  double res;
  int type_em=-1;
  Energy pt_last;
  for(int j=0;j<3;j++) {  
    pt_last=maxp;
    double a = _alphaS->overestimateValue()*_prefactor[j]*(maxyj-minyj)/(_power-1.);
    do {
      pt=GeV/pow(pow(GeV/pt_last,_power-1)-log(UseRandom::rnd())/a,1./(_power-1.));
      yj=UseRandom::rnd()*(maxyj-minyj)+ minyj;
      res=getResult(j,pt,yj);
      wgt = res/(_prefactor[j]*pow(GeV/pt,_power));
      reject = UseRandom::rnd()>wgt;
      pt_last=pt;
      //no emission event if p goes past p min - basically set to outside
      //of the histogram bounds (hopefully hist object just ignores it)
      if(pt<minp){
	pt=0.0*GeV;
	reject = false;
      }
      if(wgt>1.0)cerr<< "PROBLEM!!!!"<<endl;
    }
    while(reject);
  
    if(pt>winning_pt){
      type_em = j;
      _pt=pt;
      _yj=yj;
      winning_pt=pt;
    }
  }

  if(winning_pt/GeV<minp/GeV){ //was this an (overall) no emission event?
    _pt=0.0*GeV;
    _yj=-10;
    _yb=-10;
    type_em = 3;
  }
  (*_htype)+=double(type_em)+0.5;
  _count[type_em]++;
  return 0;
}

void DrellYanHardGenerator::dofinish() {
  //this is called at the end of run - wite histograms to hist.top
  HardestEmissionGenerator::dofinish();

  ofstream hist_out("hist3.top");
  ofstream scatxy("scatxy.dat");
  ofstream scatyjp("scatyjp.dat");

  _hyb->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yb"), string(), string(), string(), string("yb"), string());
  _hyj->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yj"), string(), string(), string(), string("yj"), string());
  _hplow->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());
  _hphigh->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());


  _htype->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("Emission type"), string(), string(), string(), string("Emission Type"), string());

  cerr<<endl<<"max(1,2,3)= "<<"("<<_max[0]<<", "<<_max[1]<<" , "<<_max[2]<<")"<<endl;

  cerr<<endl<< "no of events of each type (qqbar2Zg,qg2Zq,qbarg2Zqbar,no emiss)= ("<<_count[0]<<", "<<_count[1]<<", "<<_count[2]<<" ,"<<_count[3]<<")"<<endl;
}

void DrellYanHardGenerator::doinitrun() {
  _power=2.0;
  _max[0]=0.0;
  _max[1]=0.0;
  _max[2]=0.0;
  _prefactor[0]=5.0;
  _prefactor[1]=3.0;
  _prefactor[2]=3.0;
  _count[0]=0;
  _count[1]=0;
  _count[2]=0;
  _count[3]=0;
  _hyj= new_ptr(Histogram(-8.0,8.0,100));
  _hplow= new_ptr(Histogram(0.0,5.0,100));
  _hphigh= new_ptr(Histogram(0.0,200.0,200));
  _hyb= new_ptr(Histogram(-6.0,6.0,100));
  _htype = new_ptr(Histogram(0.0,4.0,4));
  _weighta= new_ptr(Histogram(0.0,100.0,100));
  _weightb= new_ptr(Histogram(0.0,100.0,100));
  _weightc= new_ptr(Histogram(0.0,100.0,100));
  //set up histograms in here- this is called at start of r
  HardestEmissionGenerator::doinitrun();
}



