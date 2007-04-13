// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorBosonQQbarHardGenerator class.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "VectorBosonQQbarHardGenerator.h"
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

void VectorBosonQQbarHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _prefactor << _power;
}

void VectorBosonQQbarHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _prefactor >> _power;
}

ClassDescription<VectorBosonQQbarHardGenerator> VectorBosonQQbarHardGenerator::initVectorBosonQQbarHardGenerator;
// Definition of the static class description member.

void VectorBosonQQbarHardGenerator::Init() {

  static ClassDocumentation<VectorBosonQQbarHardGenerator> documentation
    ("There is no documentation for the VectorBosonQQbarHardGenerator class");

  static Reference<VectorBosonQQbarHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VectorBosonQQbarHardGenerator::_alphaS, false, false, true, false, false);
}


NasonTreePtr VectorBosonQQbarHardGenerator::generateHardest(ShowerTreePtr tree) {
  // cerr<<"generating hardest emission"<<endl;
  int ix=0;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt){
    _quark[ix]=cjt->first->copy()->momentum();
    //swap so quark momentum is first
    if(ix==1&&cjt->first->progenitor()->id()>0) 
      swap(_quark[0],_quark[1]);
    ix++;					     
  }
  // cerr<<"before event"<<endl;

  getEvent();
  // cerr<<"after event"<<endl;
  double thrust;
  if(_qemission) thrust=_x1;
  else thrust=_x2;
  Energy mass;
  if((_pa)(3)>(_pb)(3) && (_pa)(3)>_g(3)) mass = ((_pb)+_g).mag();
  else if((_pb)(3)>(_pa)(3) && (_pb)(3)>_g(3)) mass = ((_pa)+_g).mag();
  else mass = ((_pb)+(_pa)).mag();
  // cerr<<"thrust = "<<thrust<<endl;
  //cerr<<"mass = "<<mass<<endl;
  
  _ptplot.push_back(_pt/GeV);
  _yplot.push_back(_y);
  _x1plot.push_back(_x1);
  _x2plot.push_back(_x2);
  (*_hthrust)+=1.-thrust;
  (*_hthrustlow)+=1.-thrust;
  (*_hmass)+=mass/GeV;
  (*_hy)+=_y;
  (*_hplow)+=_pt/GeV;
  (*_hphigh)+=_pt/GeV;
  //  cerr<<"finished generating hardest"<<endl;
}

bool VectorBosonQQbarHardGenerator::canHandle(ShowerTreePtr tree) {
  // two incoming particles
  // cerr<<tree->incomingLines()<<endl;
  // cerr<<"starting can handle"<<endl;
  if(tree->incomingLines().size()!=1){
    // cerr<<endl<<"not correct incoming lines = "<<tree->incomingLines().size()<<endl;
    // cerr<<"first id is "<<tree->incomingLines().begin()->first->progenitor()->id()<<endl;
    // cerr<<"first id of outgoing is "<<tree->outgoingLines().begin()->first->progenitor()->id()<<endl;
    return false;    
  }
  //  else cerr<<"passed incoming no"<<endl;

  if((tree->incomingLines().begin()->first->id()==22)&&(tree->incomingLines().begin()->first->progenitor()->id()==23)){ 
    //  cerr<<"2"<<endl;
    //   cerr<<"id is: "<<tree->incomingLines().begin()->first->progenitor()->id()<<endl;
    return false;
  }
  // else cerr<<"passes incoming no"<<endl;

  // two outgoing particles
  if((tree->outgoingLines().size()!=2)){
        
    //   cerr<<endl<<"first id is "<<tree->incomingLines().begin()->first->progenitor()->id()<<endl;
    //   cerr<<"1- out going lines = "<<tree->outgoingLines().size()<<endl; 
    // cerr<<"first id of outgoing is "<<tree->outgoingLines().begin()->first->progenitor()->id()<<endl;
    return false;
  }
  //  else cerr<<"passed outgoing"<<endl;
  // find the outgoing particles

  unsigned int ix(0);
  ShowerParticlePtr part[2];
  ix=0; 
  // cerr<<"just testing the type of outgoing"<<endl;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    {part[ix]=cjt->first->progenitor();++ix;}
  // outgoing particles check q qbar
  if(!(part[0]->id()>0&&part[0]->id()<6&&part[1]->id()<0&&part[1]->id()>-6)&&
     !(part[1]->id()>0&&part[1]->id()<6&&part[0]->id()<0&&part[0]->id()>-6)){
    //   cerr<<"4"<<endl;
    return false;
  }
  // cerr<<"can handle event"<<endl;
  return true;
}

void VectorBosonQQbarHardGenerator::doinit() throw(InitException) {
  HardestEmissionGenerator::doinit();
  // cerr<<"doinit fine"<<endl;
}

void VectorBosonQQbarHardGenerator::dofinish() {
  //this is called at the end of run - wite histograms to hist.top
  HardestEmissionGenerator::dofinish();

  ofstream hist_out("hist3.top");
  ofstream scatxy("scatxy.dat");
  ofstream scatyp("scatyp.dat");

  _hy->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist y"), string(), string(), string(), string("y"), string());
  _hplow->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());
  _hphigh->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());
  _hthrust->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Thrust"), string(), string(), string(), string("1-T"), string());
  _hthrustlow->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Thrust low"), string(), string(), string(), string("1-T"), string());
  _hmass->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Mass/GeV"), string(), string(), string(), string("Mass/GeV"), string());


  for(unsigned int i =0; i<_yplot.size();i++){
    scatxy<<_x1plot[i]<<" "<<_x2plot[i]<<endl;
    scatyp<<_yplot[i]<<" "<<_ptplot[i]<<endl;
  }
}

void VectorBosonQQbarHardGenerator::doinitrun() {
  // cerr<<"do init run"<<endl;
  _s=sqr(generator()->maximumCMEnergy());
  cerr<<"s = "<<_s/pow(GeV,2.)<<endl;
  _eventFrame = getTransf();
  // cerr<<"2"<<endl;
  _power=1.0; //remember 1 is a special case
  // cerr<<"3"<<endl;
  _prefactor=1.0;
  //cerr<<"4"<<endl;
  // getMax(10000);
    cerr<<"max is:"<<_prefactor<<endl;
  _hmass = new_ptr(Histogram(0.,60.,100));
  //cerr<<"5"<<endl;
  _hthrust = new_ptr(Histogram(0.,0.5,100));
  _hthrustlow= new_ptr(Histogram(0.0,0.01,100));
  //cerr<<"6"<<endl;
  _hy= new_ptr(Histogram(-8.0,8.0,100));
  //cerr<<"7"<<endl;
  _hplow= new_ptr(Histogram(0.0,5.,100));
  //cerr<<"8"<<endl;
  _hphigh= new_ptr(Histogram(0.0,100.0,100));
  //set up histograms in here- this is called at start of r
  // cerr<<"9"<<endl;
  HardestEmissionGenerator::doinitrun();
  // cerr<<"finish doinitrun "<<endl;
}

//private functions-internal workings




Lorentz5Momentum VectorBosonQQbarHardGenerator::getEvent(){
  
  Energy minp = 0.1*GeV;  
  Energy maxp = sqrt(0.5)*generator()->maximumCMEnergy();
  double miny = -8.0;
  double maxy = 8.0;
  double wgt;
  bool reject;
  Energy last_pt;
  last_pt=maxp;
  
  do {
    //  _pt=GeV/pow(
    //	       pow(GeV/last_pt,_power-1.)-log(UseRandom::rnd())
    //       *(_power-1.)/_prefactor/(maxy-miny),
    //       1./(_power-1.));

      _pt=last_pt*pow(UseRandom::rnd(),
    		    1./(maxy-miny)/_prefactor);//veto for c/pt case

     _y=UseRandom::rnd()*(maxy-miny)+ miny;

     _x1 = 1.-_pt/sqrt(_s)*exp(-_y);
     _x2 = 1.-_pt/sqrt(_s)*exp(_y);

     wgt = getResult() / ( _prefactor*pow(GeV/_pt,_power));
     reject = UseRandom::rnd()>wgt || ! inRange();
     
     last_pt=_pt;
     if( inRange() ){
       //  last_pt=_pt;
       if ( wgt>1.0 ){ 
	 cerr << "PROBLEM!!!!"<< endl;
	 cerr<< " res = "<< getResult() << "overfn = " << 
	   _prefactor*pow(GeV/_pt,_power) << endl;
       }
     }
     //no emission event if p goes past p min - basically set to outside
     //of the histogram bounds (hopefully hist object just ignores it)
     if(_pt<minp){
       _pt = 0.0 * GeV;//no emission event
       _y = -10;
       reject = false;
     }
  }while ( reject );

  //generate herwig variables (need to choose 1->2 splitting type)
  if(_x1>_x2){
    _qemission = true;
    _z=(_x1+_x2-1.0)/_x2;
    _ktild=_x2*_x2*(1.0-_x2)/((_x1+_x2-1.0)*(1.0-_x1));
  }
  else{
    _qemission = false;
    _z=(_x2+_x1-1.0)/_x1;
    _ktild=_x1*_x1*(1.0-_x1)/((_x2+_x1-1.0)*(1.0-_x2));
  }

  _k = sqrt(_z*_z*(1.0-_z)*(1.0-_z)*_ktild);

  //construct vectors in com z frame
  constructVectors();
  azimuthal();
 
  //boost constructed vectors into the event frame
  (_pa)=_eventFrame*(_pa);
  (_pb)=_eventFrame*(_pb);
  return _g=_eventFrame*_g;
  }

 double VectorBosonQQbarHardGenerator::getResult() {
   //   cerr<<"_x1,_x2,_pt,y= "<<_x1<<", "<<_x2<<", "<<_pt<<", "<<_y<<endl;
   //   double res = 4./3./pi*2./_pt*
   //  (1.-2.*_pt*cosh(_y)/sqrt(_s)+sqr(_pt)*cosh(2.*_y)/_s)*GeV;
   double res=4./3./pi*_pt/_s*(sqr(_x1)+sqr(_x2))/(1.-_x1)/(1.-_x2)*GeV;
   //cerr<<res<<endl;
   res*= _alphaS->value(sqr(_pt));
   return res;
 }

double VectorBosonQQbarHardGenerator::getMax(int num){
  
  double res;   
  Energy minp = 0.1*GeV;  
  Energy maxp = sqrt(0.5)*generator()->maximumCMEnergy();
  double miny = -8.;
  double maxy = 8.;
  _max=0.;

  for(int i =0;i<num;i++){
    do{
      _pt = UseRandom::rnd()*(maxp-minp)+ minp;
      _y = UseRandom::rnd()*(maxy-miny)+ miny;
      _x1 = 1.-_pt/sqrt(_s)*exp(-_y);
      _x2 = 1.-_pt/sqrt(_s)*exp(_y);
      res = getResult();
      //    cerr<<"res= "<<res<<endl;;
    }while(!inRange());
    if (res*pow((_pt/GeV),_power)>_max) _max = res*pow((_pt/GeV),_power);
  }

  cerr<<"max is: "<<_max<<endl;
  _prefactor=_max;
  return _max;
}



//momentum construction
LorentzRotation VectorBosonQQbarHardGenerator::getTransf(){
 
  LorentzRotation transf((_quark[0]+_quark[1]).findBoostToCM());

  Lorentz5Momentum q1=transf*(_pa);
  
  if(q1(0)==0.0) transf.rotateZ(-pi/2.0);
  else transf.rotateZ(-atan(q1(1)/q1(0)));

  if(q1(2)==0.0) transf.rotateY(pi/2.0);
  else transf.rotateY(atan(sqrt(q1(0)*q1(0)+q1(1)*q1(1))/q1(2)));

  Lorentz5Momentum q2=transf*(_pa);

  if(q2(2)<0.0)transf.rotateY(pi);
  
  transf.invert();

  return transf;
}


void VectorBosonQQbarHardGenerator::azimuthal() {
   if(UseRandom::rnd()<_x1*_x1/(_x1*_x1+_x2*_x2)){
     _r.setRotate(UseRandom::rnd()*2.0*pi, (_pa).vect() );
     (_pb)=_r *(_pb);
   }
   else{
     _r.setRotate( UseRandom::rnd()*2.0*pi, (_pb).vect() );
     _pa=_r *(_pa);
   }
    _g=_r*_g;
    return;
}


void VectorBosonQQbarHardGenerator::constructVectors(){

  _phi=UseRandom::rnd()*pi;
  if(_qemission){
   (_pa).setT(sqrt(_s)*(_z+_k*_k/_z)/2.0);
   (_pa).setX(sqrt(_s)*_k*cos(_phi));
   (_pa).setY(sqrt(_s)*_k*sin(_phi));
   (_pa).setZ(sqrt(_s)*(_z-_k*_k/_z)/2.0);

   (_pb).setT(sqrt(_s)*(1.0-_k*_k/_z/(1.0-_z))/2.0);
   (_pb).setX(0.0);
   (_pb).setY(0.0);
   (_pb).setZ(sqrt(_s)*(_k*_k/_z/(1-_z)-1.0)/2.0);
    
   _g.setT(sqrt(_s)*(1.0-_z+_k*_k/(1.0-_z))/2.0);
   _g.setX(sqrt(_s)*-_k*cos(_phi));
   _g.setY(sqrt(_s)*-_k*sin(_phi));
   _g.setZ(sqrt(_s)*(1.0-_z-_k*_k/(1.0-_z))/2.0);
  }
  else{
   (_pa).setT(sqrt(_s)*(1.0-_k*_k/_z/(1.0-_z))/2.0);
   (_pa).setX(0.0);
   (_pa).setY(0.0);
   (_pa).setZ(sqrt(_s)*(1.0-_k*_k/_z/(1.0-_z))/2.0);

   (_pb).setT(sqrt(_s)*(_z+_k*_k/_z)/2.0);
   (_pb).setX(sqrt(_s)*_k*cos(_phi));
   (_pb).setY(sqrt(_s)*_k*sin(_phi));
   (_pb).setZ(sqrt(_s)*(-_z+_k*_k/_z)/2.0);
    
   _g.setT(sqrt(_s)*((1.0-_z)+_k*_k/(1.0-_z))/2.0);
   _g.setX(sqrt(_s)*-_k*cos(_phi));
   _g.setY(sqrt(_s)*-_k*sin(_phi));
   _g.setZ(sqrt(_s)*(-(1.0-_z)+_k*_k/(1.0-_z))/2.0);
  }
  return;  
}




