// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanMECorrection class.
//

#include "DrellYanMECorrection.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DrellYanMECorrection.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DrellYanMECorrection::~DrellYanMECorrection() {}

void DrellYanMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _phasespaceopt << _channelwgt << _channelweights;
}

void DrellYanMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _phasespaceopt >> _channelwgt >> _channelweights;
}

ClassDescription<DrellYanMECorrection> DrellYanMECorrection::initDrellYanMECorrection;
// Definition of the static class description member.

void DrellYanMECorrection::Init() {

  static ClassDocumentation<DrellYanMECorrection> documentation
    ("The DrellYanMECorrection class implements the soft and hard"
     " matrix element correction for the Drell-Yan process. This is"
     " a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.");
  static Switch<DrellYanMECorrection,unsigned int> interfacePhaseSpaceOption
    ("PhaseSpaceOption",
     "Option for the generation of the phase-space",
     &DrellYanMECorrection::_phasespaceopt, 0, false, false);
  static SwitchOption interfacePhaseSpaceOptionMandelstam
    (interfacePhaseSpaceOption,
     "Mandelstam",
     "Generate the phase space in the Mandelstam variables",
     0);
  static SwitchOption interfacePhaseSpaceOptionRapidityAndPt
    (interfacePhaseSpaceOption,
     "RapidityAndPt",
     "Generate the phase space in rapidity and transverse momentum",
     1);

  static Parameter<DrellYanMECorrection,double> interfaceChannelWeight
    ("ChannelWeight",
     "The relative weights of the q qbar abd q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &DrellYanMECorrection::_channelwgt, 0.1, 0.01, 100.,
     false, false, Interface::limited);

}

bool DrellYanMECorrection::canHandle(ShowerTreePtr tree)
{
   // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be a quark and an antiquark
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
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
  for(cit=tree->outgoingLines().begin();cit!=tree->outgoingLines().end();++cit)
    {part[ix]=cit->first->progenitor();++ix;}
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
  cerr << "testing apply drell yan " << endl;
  return true;
}

void DrellYanMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree)
{
  cerr << "testing start of hard " << endl;
  // get the quark,antiquark and the gauge boson
  // get the quarks
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  ShowerParticleVector incoming;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit)
    incoming.push_back(cit->first->progenitor());
  // get the gauge bosons
  PPtr boson;
  if(tree->outgoingLines().size()==1)
    {boson=tree->outgoingLines().begin()->first->copy();}
  else
    {boson=tree->outgoingLines().begin()->first->copy()->parents()[0];}
  // ensure that the quark is first
  bool quarkfirst(true);
  if(incoming[0]->id()<incoming[1]->id())
    {
      quarkfirst=false;
      swap(incoming[0],incoming[1]);
    }
  unsigned int iemit,itype;
  vector<Lorentz5Momentum> pnew;
  if(!applyHard(incoming,boson,iemit,itype,pnew));
}

bool DrellYanMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
						 ShowerParticlePtr parent,Branching br)
{
//   //cerr << "testing in top decay soft correction " 
//   //     << *initial->progenitor() << "\n" << *parent << "\n";
//   // check if we need to apply the full correction
//   unsigned int id[2]={abs(initial->progenitor()->id()),abs(parent->id())};
//   // the initial-state correction
//   if(id[0]==ParticleID::t&&id[1]==ParticleID::t)
//     {
//       Energy pt=br.kinematics->pT();
//       // check if hardest so far
//       // if not just need to remove effect of enhancement
//       bool veto(false);
//       // if not hardest so far
//       if(pt<initial->pT())
// 	veto=!UseRandom::rndbool(1./_initialenhance);
//       // if hardest so far do calculation
//       else
// 	{
// 	  // values of kappa and z
// 	  double z(br.kinematics->z()),kappa(sqr(br.kinematics->qtilde()/_mt));
// 	  // parameters for the translation
// 	  double w(1.-(1.-z)*(kappa-1.)),u(1.+_a-_c-(1.-z)*kappa),v(sqr(u)-4.*_a*w*z);
// 	  // veto if outside phase space
// 	  if(v<0.) 
// 	    veto=true;
// 	  // otherwise calculate the weight
// 	  else
// 	    {
// 	      v = sqrt(v);
// 	      double xa((0.5*(u+v)/w+0.5*(u-v)/z)),xg((1.-z)*kappa);
// 	      double f(me(xa,xg)),
// 		J(0.5*(u+v)/sqr(w)-0.5*(u-v)/sqr(z)+_a*sqr(w-z)/(v*w*z));
// 	      double wgt(f*J*2./kappa/(1.+sqr(z)-2.*z/kappa)/_initialenhance);
// 	      if(wgt>1.)
// 		{
// 		  cerr << "testing violates max initial " << xg 
// 		       << " " << xa << " " << wgt << " " << _initialenhance << endl;
// 		  wgt=1.;
// 		}
// 	      // compute veto from weight
// 	      veto = !UseRandom::rndbool(wgt);
// 	    }
// 	  // if not vetoed reset max
// 	  if(!veto) initial->pT(pt);
// 	}
//       // if vetoing reset the scale
//       if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->qtilde());
//       // return the veto
//       return veto;
//     }
//   // final-state correction
//   else if(id[0]==ParticleID::b&&id[1]==ParticleID::b)
//     {
//       Energy pt=br.kinematics->pT();
//       // check if hardest so far
//       // if not just need to remove effect of enhancement
//       bool veto(false);
//       // if not hardest so far
//       if(pt<initial->pT()) return !UseRandom::rndbool(1./_finalenhance);
//       // if hardest so far do calculation
//       // values of kappa and z
//       double z(br.kinematics->z()),kappa(sqr(br.kinematics->qtilde()/_mt));
//       // momentum fractions
//       double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),
// 	xg((2.-xa)*(1.-r)-(z-r)*(sqr(xa)-4.*_a));
//       double root(sqr(xa)-4.*_a),xfact(z*kappa/2./(z*(1.-z)*kappa+_c)*(2.-xa-root)+root);
//       // calculate the full result
//       double f(me(xa,xg));
//       // jacobian
//       double J(z*root);
//       double wgt(f*J*2.*kappa/(1.+sqr(z)-2.*_c/kappa/z)/sqr(xfact)/_finalenhance);
//       if(wgt>1.)
// 	{
// 	  cerr << "testing violates max final " << xg 
// 	       << " " << xa << " " << wgt << endl;
// 	  wgt=1.;
// 	}
//       // compute veto from weight
//       veto = !UseRandom::rndbool(wgt);
//       // if not vetoed reset max
//       if(!veto) initial->pT(pt);
//       // if vetoing reset the scale
//       if(veto) parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->qtilde());
//       // return the veto
//       return veto;
//     }
//   // otherwise don't veto
//   else return !UseRandom::rndbool(1./_finalenhance);
  return false;
}

bool DrellYanMECorrection::applyHard(ShowerParticleVector quarks, PPtr boson,
				     unsigned int & iemit,unsigned int & itype,
				     vector<Lorentz5Momentum> & pnew)
{
  // calculate the limits on s
  Energy mb(boson->mass());
  Energy2 smin=sqr(mb),mb2(smin);
  Energy2 s=
    (CurrentGenerator::current().currentEvent()->incoming().first->momentum()+
     CurrentGenerator::current().currentEvent()->incoming().second->momentum()).m2();
  Energy2 smax(s);
  // calculate the rapidity of the boson
  double yB=0.5*log((boson->momentum().e()+boson->momentum().pz())/
		    (boson->momentum().e()-boson->momentum().pz()));
  // if no phase-space return
  if(smax<smin) return false;
  // get the evolution scales (this needs improving)
  double kappa[2]={1.,1.};
  // get the momentum fractions for the leading order process
  // and the values of the PDF's
  double x[2],fx[2];
  tcPDFPtr pdf[2];
  for(unsigned int ix=0;ix<quarks.size();++ix)
    {
      x[ix]=quarks[ix]->x();
      Ptr<BeamParticleData>::const_pointer beam=
	dynamic_ptr_cast<Ptr<BeamParticleData>::const_pointer>
	(quarks[ix]->parents()[0]->dataPtr());
      assert(beam);
      pdf[ix]=beam->pdf();
      assert(pdf[ix]);
      fx[ix]=pdf[ix]->xfx(quarks[ix]->parents()[0]->dataPtr(),quarks[ix]->dataPtr(),
			  mb2,x[ix]);
    }
  // select the type of process and generate the kinematics
  double rn(UseRandom::rnd());
  Energy shat,uhat,that;
  double weight(0.),xnew[2];
  cerr << "testing B " << _channelweights.size() << endl;
  // q qbar -> g V
  if(rn<_channelweights[0])
    {
      double jacobian(1.);
      // generate in s and t
      if(_phasespaceopt==0)
	{
	  cerr << "testing in A " << endl;
	  // generate the value of s according to 1/s^2
	  shat = smax*smin/(smin+UseRandom::rnd()*(smax-smin));
	  cerr << "testing shat " << smin << " " << smax << " " << shat << endl;
	  double sbar=shat/mb2;
	  jacobian = sqr(shat)*(1./smin-1./smax); 
	  // calculate the limits on that and uhat
	  Energy tmax=mb2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
	  Energy tmin=shat*(1.-sbar)/(kappa[1]+sbar);
	  Energy umax(mb2-shat-tmin),umin(mb2-shat-tmax);
	  // check inside phase space
	  if(tmax<tmin||umax<umin) return false;
	  cerr << "testing t limits " << tmax << " " << tmin << endl;
	  cerr << "testing u limits " << umax << " " << umin << endl;
	  // generate t and u according to 1/t+1/u
	  // generate in 1/t
	  if(UseRandom::rndbool(0.5))
	    {
	      that=tmax*pow(tmin/tmax,UseRandom::rnd());
	      uhat=mb2-shat-that;
	      jacobian *=log(tmin/tmax);
	    }
	  // generate in 1/u
	  else
	    {
	      uhat=umax*pow(umin/umax,UseRandom::rnd());
	      that=mb2-shat-uhat;
	      jacobian *=log(umin/umax);
	    }
	  jacobian *=2.*uhat*that/(shat-mb2);
	}
      // generate in y and pt
      else
	{
	}
      // new scale (this is mt^2=pt^2+mb^2)
      Energy2 scale(uhat*that/shat+mb2);
      // the PDF's with the emitted gluon
      double xnew[2],fxnew[2];
      xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(mb2-uhat)/(mb2-that));
      xnew[1]=shat/(s*xnew[0]);
      for(unsigned int ix=0;ix<2;++ix)
	{fxnew[ix]=pdf[ix]->xfx(quarks[ix]->parents()[0]->dataPtr(),
				quarks[ix]->dataPtr(),mb2,xnew[ix]);}
      // jacobian and me parts of the weight
      weight=jacobian*(sqr(mb2-that)+sqr(mb2-uhat))/(sqr(shat)*that*uhat);
      // pdf part of the weight
      weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
      // finally coupling, colour factor and different channel pieces
      weight *= 2./3./pi/_channelweights[0]*coupling()->value(scale);
      cerr << "testing the weight " << weight << " " << shat << " " << that << endl;
      // select the emiting particle
      iemit=1;
      if(UseRandom::rnd()<sqr(mb2-uhat)/(sqr(mb2-uhat)+sqr(mb2-that))) iemit=2;
      itype=0;
    }
  // incoming gluon
  else
    {
      double jacobian;
      // generate in s and t 
      if(_phasespaceopt==0)
 	{
 	  // generate the value of s according to 1/s^2
 	  shat = smax*smin/(smin+UseRandom::rnd()*(smax-smin));
 	  jacobian = sqr(shat)*(1./smin-1./smax);
 	  double sbar=shat/mb2;
	  // generate t 
	  // calculate limits
 	  Energy tmax=mb2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
	  Energy tmin=shat*(1.-sbar)/(kappa[1]+sbar);
	  if(rn>_channelweights[1])
	    {
	      swap(tmax,tmin);
	      tmax=mb2-shat-tmax;
	      tmin=mb2-shat-tmin;
	    }
 	  // check inside phase space
 	  if(tmax<tmin) return false;
	  that=tmax*pow(tmin/tmax,UseRandom::rnd());
	  uhat=mb2-shat-that;
	  jacobian *=that*log(tmax/tmin);
	}
      else
	{
	}
      // new scale (this is mt^2=pt^2+mb^2)
      Energy2 scale(uhat*that/shat+mb2);
      // g qbar -> qbar V 
      double fxnew[2];
      if(rn<_channelweights[1])
	{
	  itype=2;
	  xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(mb2-uhat)/(mb2-that));
	  xnew[1]=shat/(s*xnew[0]);
	  fxnew[0]=pdf[0]->xfx(quarks[0]->parents()[0]->dataPtr(),
			       getParticleData(ParticleID::g),mb2,xnew[0]);
	  fxnew[1]=pdf[1]->xfx(quarks[1]->parents()[0]->dataPtr(),
			       quarks[1]->dataPtr(),mb2,xnew[1]);
	  jacobian/=(_channelweights[1]-_channelweights[0]);
	}
      // q g -> q V 
      else
	{
	  itype=1;
	  xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(mb2-that)/(mb2-uhat));
	  xnew[1]=shat/(s*xnew[0]);
	  fxnew[0]=pdf[0]->xfx(quarks[0]->parents()[0]->dataPtr(),
			       quarks[0]->dataPtr(),mb2,xnew[0]);
	  fxnew[1]=pdf[1]->xfx(quarks[1]->parents()[0]->dataPtr(),
			       getParticleData(ParticleID::g),mb2,xnew[1]);
	  jacobian/=(_channelweights[2]-_channelweights[1]);
	}
      // jacobian and me parts of the weight
      weight=-jacobian*(sqr(mb2-that)+sqr(mb2-shat))/(sqr(shat)*shat*that);
      // pdf part of the weight
      weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
      // finally coupling, colour factor and different channel pieces
      weight *= 0.25/pi*coupling()->value(scale);
      cerr << "testing the weight " << weight << " " << shat << " " << that << endl;
      // select the emiting particle
      iemit=1;
      if(UseRandom::rnd()<sqr(mb2-that)/(sqr(mb2-that)+sqr(mb2-shat))) iemit=2;
    }
  // if me correction should be applied
  if(weight>1.)
    {
      cerr << "testinmg weight greater than 1 " << weight;
      weight=1.;
    }
  if(UseRandom::rnd()>weight) return false;
  // construct the momenta in the rest frame of the boson
  Lorentz5Momentum pb(0.,0.,0.,mb,mb),pspect,pg,pemit;
  double cos3;
  if(itype==0)
    {
      pg     = Lorentz5Momentum(0.,0.,0.,0.5*(shat-mb2)/mb,0.);
      Energy2 tp(that),up(uhat);
      double zsign(-1.);
      if(iemit==2)
	{
	  tp=uhat;
	  up=that;
	  zsign=1.;
	}
      pspect = Lorentz5Momentum(0.,0.,zsign*0.5*(mb2-tp)/mb,0.5*(mb2-tp)/mb,0.);
      Energy eemit=0.5*(mb2-up)/mb;
      cos3 = 0.5/pspect.pz()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
  else
    {
      pg=Lorentz5Momentum(0.,0.,0.,0.5*(mb2-uhat)/mb,0.);
      double zsign(1.);
      if(iemit==1)
	{
	  if(itype==1) zsign=-1.;
	  pspect=Lorentz5Momentum(0.,0.,0.5*zsign*(shat-mb2)/mb,0.5*(shat-mb2)/mb);
	  Energy eemit=0.5*(mb2-that)/mb;
	  cos3 = 0.5/pspect.pz()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
	}
      else
	{
	  if(itype==2) zsign=-1.;
	  pspect=Lorentz5Momentum(0.,0.,0.5*zsign*(mb2-that)/mb,0.5*(mb2-that)/mb);
	  Energy eemit=0.5*(shat-mb2)/mb;
	  cos3 = 0.5/pspect.pz()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
	}
    }
  // rotate the gluon
  double sin3(sqrt(1.-sqr(cos3))),phi(2.*pi*UseRandom::rnd());
  pg.setPx(pg.e()*sin3*cos(phi));
  pg.setPy(pg.e()*sin3*sin(phi));
  pg.setPz(pg.e()*cos3);
  if(itype==0)
    {pemit=pb+pg-pspect;}
  else
    {
      if(iemit==1)
	pemit=pb+pspect-pg;
      else
	pemit=pspect+pg-pb;
    }
  pemit.rescaleMass();
  if(itype==0)
    {
      cerr << "testing momenta " 
	   << pb     << " " << pb.m()     << "\n" 
	   << pspect << " " << pspect.m() << "\n" 
	   << pg     << " " << pg.m()     << "\n"
	   << pemit  << " " << pemit.m()  << "\n";
      cerr << "testing cons " << pspect+pemit-pb-pg << endl;
    }
  else
    {
      cerr << "testing momenta " 
	   << pb     << " " << pb.m()     << "\n" 
	   << pspect << " " << pspect.m() << "\n" 
	   << pg     << " " << pg.m()     << "\n"
	   << pemit  << " " << pemit.m()  << "\n";
      cerr << "testing cons " << pspect+pg-pb-pemit << endl;
    }
  cerr << "testing weight !!! " << weight << endl;


  // find the new CMF
  Lorentz5Momentum pp[2];
  if(itype==0)
    {
      if(iemit==1)
	{
	  pp[0]=pemit;
	  pp[1]=pspect;
	}
      else
	{
	  pp[0]=pspect;
	  pp[1]=pemit;
	}
    }
  else if(itype==1)
    {
      pp[1]=pg;
      if(iemit==1) pp[0]=pemit;
      else         pp[0]=pspect;
    }
  else
    {
      pp[0]=pg;
      if(iemit==1) pp[1]=pemit;
      else         pp[1]=pspect;
    }
  Lorentz5Momentum pz=pp[0];
  pp[0]/=xnew[0];
  pp[1]/=xnew[1];
  Lorentz5Momentum plab=pp[0]+pp[1];
  plab.rescaleMass();
  // construct the boost to rest frame of plab
  LorentzRotation boost(plab.findBoostToCM());
  pz.transform(boost); 
  // rotate so emitting particle along z axis
  LorentzRotation rot;
  // rotate so in x-z plane
  rot.setRotateZ(-atan2(pz.y(),pz.x()));
  rot.rotateY(-acos(pz.z()/pz.vect().mag()));
  rot.rotateZ(atan2(pz.y(),pz.x()));
  pz.transform(rot);
  cerr << "testing pz " << pz << endl;

  exit(0);
  return true;
}
