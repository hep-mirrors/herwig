// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleLHCAnalysis class.
//

#include "SimpleILCAnalysis.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SimpleILCAnalysis.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SimpleILCAnalysis::SimpleILCAnalysis() :
  _acostem(0.,3.2,32),
  _acostee(0.,3.2,32),
  _acostebe(0.,3.2,32),
  _Mt(120.,200.,100),
  _Mtb(120.,200.,100),
  _rapb(-5.,5.,50),
  _rapbb(-5.,5.,50),
  _pzb(-200.,200.,80),
  _pzbb(-200.,200.,80),
  _ptb(0.,200.,40),
  _ptbb(0.,200.,40),
  _pbe(0.,200.,40),
  _pbbe(0.,200.,40),
  _pzbr(0.,1.,50),
  _pzbbr(0.,1.,50),
  _ptbr(0.,1.,50),
  _ptbbr(0.,1.,50)

{}

   
void SimpleILCAnalysis::analyze(tEventPtr event, long, int, int) {
   eventweight_=event->weight();
   cout << "eventweight" << "                " << eventweight_ << endl;
   // Rotate to CMS, extract final state particles and call analyze(particles).
  
  Lorentz5Momentum pz;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator send=event->primaryCollision()->steps().end();
  
  ParticleVector b;
  b.clear();
  ParticleVector lept;
  lept.clear();

  Lorentz5Momentum pb(0.0,0.0,0.0,0.0,0.0);
  Lorentz5Momentum pbb(0.0,0.0,0.0,0.0,0.0);
 
  for( ;sit!=send;sit=send)
   
  for(sit=send-2;sit!=send-1;++sit)          {
   ParticleSet partt=(**sit).all();
   ParticleSet::const_iterator iterr=partt.begin();
   ParticleSet::const_iterator endd =partt.end();
   for( ;iterr!=endd;++iterr) {
	if (fabs((**iterr).id())==5) {
	 
      b.push_back((*iterr));
    }
   }
  ParticleVector::const_iterator iters=b.begin();
  ParticleVector::const_iterator ends =b.end();
  for( ;iters!=ends;++iters) {
      if ((**iters).id()==5) {
      pb=(*iters)->momentum();
          } else 
      if ((**iters).id()==-5) {
      pbb=(*iters)->momentum();  
   
  }
  }
  }
 
   for(sit=send-2;sit!=send-1;++sit)
   {
   ParticleSet part=(**sit).all();
   ParticleSet::const_iterator iter=part.begin();
   ParticleSet::const_iterator end =part.end();

   for( ;iter!=end;++iter) {
        if (((fabs((**iter).id())==11) || (fabs((**iter).id())==12) )&& (**iter).children().size()==0) {
          lept.push_back((*iter));
     }
   }
  

   Lorentz5Momentum pep(0.0,0.0,0.0,0.0,0.0);
   Lorentz5Momentum pem(0.0,0.0,0.0,0.0,0.0);
   Lorentz5Momentum pnup(0.0,0.0,0.0,0.0,0.0);
   Lorentz5Momentum pnum(0.0,0.0,0.0,0.0,0.0);
   ParticleVector::const_iterator iter1=lept.begin();
   ParticleVector::const_iterator end1 =lept.end();

   for( ;iter1!=end1;++iter1) {
     if((**iter1).id()==-11) {
       pep=(*iter1)->momentum();}
     else if ((**iter1).id()==11) {
       pem=(*iter1)->momentum();}
     else if ((**iter1).id()==-12) {
       pnup=(*iter1)->momentum();}
     else if ((**iter1).id()==12) {
       pnum=(*iter1)->momentum();}
    
   }
 
  Energy pbe=pb.e()/GeV;
  Energy pbbe=pbb.e()/GeV;
  Energy ptb=pb.perp()/GeV;
  Energy ptbb=pbb.perp()/GeV;
  double ptbr=ptb/pbe;
  double ptbbr=ptbb/pbbe;
 
   Energy pzb=pb.z()/GeV;
   Energy pzbb=pbb.z()/GeV;
   double pzbr=pzb/pbe;
   double pzbbr=pzbb/pbbe;
   double rapb=0.5*log((pbe+pzb)/(pbe-pzb));
   double rapbb=0.5*log((pbbe+pzbb)/(pbbe-pzbb));
   Lorentz5Momentum pt=pb+pep+pnum;
   Lorentz5Momentum ptbar=pbb+pem+pnup;
  Energy Mt=pt.mass()/GeV;
  Energy Mtb=ptbar.mass()/GeV;

 
      double costem = (pem.x()*pt.x()+pem.y()*pt.y()+pem.z()*pt.z())/(sqrt(pem.perp()*pem.perp()+pem.z()*pem.z())*sqrt(pt.perp()*pt.perp()+pt.z()*pt.z()));
       double acostem = acos(costem);
      
       double costee = (pem.x()*pep.x()+pem.y()*pep.y()+pem.z()*pep.z())/(sqrt(pem.x()*pem.x()+pem.y()*pem.y()+pem.z()*pem.z())*sqrt(pep.perp()*pep.perp()+pep.z()*pep.z()));
       double acostee = acos(costee);
     
       double costebe= pem.z()/sqrt(pem.perp()*pem.perp()+pem.z()*pem.z());
       double acostebe = acos(costebe);
       
     	_acostem.addWeighted(acostem,eventweight_);
 	_acostee.addWeighted(acostee,eventweight_);
 	_acostebe.addWeighted(acostebe,eventweight_);
        _Mt.addWeighted(Mt,eventweight_);
	_Mtb.addWeighted(Mtb,eventweight_);
	_pbe.addWeighted(pbe,eventweight_);
	_pbbe.addWeighted(pbbe,eventweight_);
	_ptb.addWeighted(ptb,eventweight_);
	_ptbb.addWeighted(ptbb,eventweight_);
	_ptbr.addWeighted(ptbr,eventweight_);
	_ptbbr.addWeighted(ptbbr,eventweight_);
	_pzbr.addWeighted(pzbr,eventweight_);
	_pzbbr.addWeighted(pzbbr,eventweight_);
       	_rapb.addWeighted(rapb,eventweight_);
       	_rapbb.addWeighted(rapbb,eventweight_);
	_pzb.addWeighted(pzb,eventweight_);
	_pzbb.addWeighted(pzbb,eventweight_);

   }

}

//LorentzRotation SimpleILCAnalysis::transform(tEventPtr) const {
//  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
//}

//void SimpleILCAnalysis::analyze(const tPVector & particles) {
//  AnalysisHandler::analyze(particles);
//}

//void SimpleILCAnalysis::analyze(tPPtr) {}

//void SimpleILCAnalysis::persistentOutput(PersistentOStream &) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
//}

//void SimpleILCAnalysis::persistentInput(PersistentIStream &, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
//}

NoPIOClassDescription<SimpleILCAnalysis> SimpleILCAnalysis::initSimpleILCAnalysis;
// Definition of the static class description member.

void SimpleILCAnalysis::Init() {

  static ClassDocumentation<SimpleILCAnalysis> documentation
    ("There is no documentation for the SimpleILCAnalysis class");

}

 void SimpleILCAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  string title;
  using namespace HistogramOptions;
 
     _Mt.topdrawOutput(outfile,Frame,"BLACK","Mass of top (e+,nu,b)/GeV");

    _Mtb.topdrawOutput(outfile,Frame,"BLACK","Mass of anti-top (e-,nubar,bbar)/GeV");

   _acostem.topdrawOutput(outfile,Frame,"BLACK","Angle between top and decay e-");
   _acostee.topdrawOutput(outfile,Frame,"BLACK","Angle between decay e+ and e-");
   _acostebe.topdrawOutput(outfile,Frame,"BLACK","Angle between decay e- and beam e-");

  _rapb.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of b quark");
   _rapbb.topdrawOutput(outfile,Frame|Ylog,"BLACK","Rapidity of anti-b quark");

   _pbe.topdrawOutput(outfile,Frame,"BLACK","Energy of b quark");
   _pbbe.topdrawOutput(outfile,Frame,"BLACK","Energy of anti-b quark");
   _ptb.topdrawOutput(outfile,Frame,"BLACK","Transverse mmt of b quark");
   _ptbb.topdrawOutput(outfile,Frame,"BLACK","Transverse mmt of anti-b quar");
   _ptbr.topdrawOutput(outfile,Frame,"BLACK","Transverse mmt/Energy ratio of b quark");
   _ptbbr.topdrawOutput(outfile,Frame,"BLACK","Transverse mmt/Energy ratio of anti-b quark");
   _pzb.topdrawOutput(outfile,Frame,"BLACK","Longitudinal mmt of b quark");
   _pzbb.topdrawOutput(outfile,Frame,"BLACK","Longitudinal mmt of anti-b quark");
   _pzbr.topdrawOutput(outfile,Frame,"BLACK","Longitudinal mmt/Energy ratio of b quark");
   _pzbbr.topdrawOutput(outfile,Frame,"BLACK","Longitudinal mmt/Energy ratio of anti-b quark");
  
}


// void SimpleILCAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
//   //  AnalysisHandler::analyze(event, ieve, loop, state);
//   // Rotate to CMS, extract final state particles and call analyze(particles).
//   // find the Z
//   Lorentz5Momentum pz;
//   StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
//   StepVector::const_iterator send=event->primaryCollision()->steps().end();
//   for(;sit!=send;++sit)
//     {
//       ParticleSet part=(**sit).all();
//       ParticleSet::const_iterator iter=part.begin();
//       ParticleSet::const_iterator end =part.end();
//       for( ;iter!=end;++iter)
// 	{
// 	  if(((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma)&&
// 	     (**iter).children().size()==2)
// 	    {
// 	      pz=getMomentum(*iter);
// 	      Energy pt = pz.perp()/GeV;
// 	      Energy mz = pz.mass()/GeV;
// 	      if(mz>20.&&mz<80.) (_ptZ[1])+=(pt);
// 	      else if (mz>80.&&mz<100.) (_ptZ[2])+=(pt);
// 	      else if (mz>100.) (_ptZ[3])+=(pt);
// 	      (_ptZ[0])+=(pt);
// 	      (_mZ)+=(mz);
// 	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
// 	      (_rapZ)+=(rap);
// 	      (_phiZ)+=pz.phi();
// 	    } else if ((**iter).id()==ParticleID::Wplus) {
// 	      pz=getMomentum(*iter);
// 	      Energy pt = pz.perp()/GeV;
// 	      Energy mz = pz.mass()/GeV;
// 	      if(mz>20.&&mz<80.) (_ptWp[1])+=(pt);
// 	      else if (mz>80.&&mz<100.) (_ptWp[2])+=(pt);
// 	      else if (mz>100.) (_ptWp[3])+=(pt);
// 	      (_ptWp[0])+=(pt);
// 	      (_mWp)+=(mz);
// 	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
// 	      (_rapWp)+=(rap);
// 	      (_phiWp)+=pz.phi();
// 	    } else if ((**iter).id()==ParticleID::Wminus) {
// 	      pz=getMomentum(*iter);
// 	      Energy pt = pz.perp()/GeV;
// 	      Energy mz = pz.mass()/GeV;
// 	      if(mz>20.&&mz<80.) (_ptWm[1])+=(pt);
// 	      else if (mz>80.&&mz<100.) (_ptWm[2])+=(pt);
// 	      else if (mz>100.) (_ptWm[3])+=(pt);
// 	      (_ptWm[0])+=(pt);
// 	      (_mWm)+=(mz);
// 	      double rap = 0.5*log((pz.e()+pz.z())/(pz.e()-pz.z()));
// 	      (_rapWm)+=(rap);
// 	      (_phiWm)+=pz.phi();
// 	    }
// 	}
//     }
// }

// LorentzRotation SimpleILCAnalysis::transform(tEventPtr event) const {
//   return LorentzRotation();
//   // Return the Rotation to the frame in which you want to perform the analysis.
// }

// void SimpleILCAnalysis::analyze(const tPVector & particles) {
//   AnalysisHandler::analyze(particles);
// }

// void SimpleILCAnalysis::analyze(tPPtr) {}

// void SimpleILCAnalysis::persistentOutput(PersistentOStream & os) const {
//   // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
// }

// void SimpleILCAnalysis::persistentInput(PersistentIStream & is, int) {
//   // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
// }

// ClassDescription<SimpleILCAnalysis> SimpleILCAnalysis::initSimpleILCAnalysis;
// // Definition of the static class description member.

// void SimpleILCAnalysis::Init() {

//   static ClassDocumentation<SimpleILCAnalysis> documentation
//     ("There is no documentation for the SimpleILCAnalysis class");

// }

