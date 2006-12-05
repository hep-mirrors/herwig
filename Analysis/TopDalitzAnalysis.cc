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
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

TopDalitzAnalysis::~TopDalitzAnalysis() {}

void TopDalitzAnalysis::analyze(tEventPtr event, long, int, int) {
  // Gets all the particles in the primaryCollision step(1)
  ParticleSet pert=event->primaryCollision()->step(1)->all();
  // Gets just the final state particles of the primaryCollision step(1)
  tPVector final=event->primaryCollision()->step(1)->getFinalState();
  ParticleSet::const_iterator pit;
  // Find the two tops from the primary collision step(1) and 
  // call topShower on each of them...
  tPVector tShower,tbarShower;
  for(pit=pert.begin();pit!=pert.end();++pit)
    {
      // All kinds of stuff is coming through here, b's, W's, t's, g's
      // c's, e's gamma's... 
      // Must be top
      if((*pit)->id()!=6) continue;
      // must have two children
      if((*pit)->children().size()!=2) continue;
      // neither should be top
      if(abs((*pit)->children()[0]->id())==6||
	 abs((*pit)->children()[1]->id())==6) continue;
      tShower=particleID(*pit,final);
      dalitz(tShower);
    }
  for(pit=pert.begin();pit!=pert.end();++pit)
    {
      // All kinds of stuff is coming through here, b's, W's, t's, g's
      // c's, e's gamma's... 
      // Must be anti-top
      if((*pit)->id()!=-6) continue;
      // must have two children
      if((*pit)->children().size()!=2) continue;
      // neither should be top
      if(abs((*pit)->children()[0]->id())==6||
	 abs((*pit)->children()[1]->id())==6) continue;
      tbarShower=particleID(*pit,final);
      dalitz(tbarShower);
    }
  Energy2 s = sqr(event->incoming().first->momentum().e()
            +     event->incoming().second->momentum().e()
                 );
  if(s!=1.296e+11) cout << "\n\ns : " << s << "\n\n";
  threeJetAnalysis(s,tShower,tbarShower);
}

LorentzRotation TopDalitzAnalysis::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void TopDalitzAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void TopDalitzAnalysis::analyze(tPPtr) {}

void TopDalitzAnalysis::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void TopDalitzAnalysis::persistentInput(PersistentIStream & , int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<TopDalitzAnalysis> TopDalitzAnalysis::initTopDalitzAnalysis;
// Definition of the static class description member.

void TopDalitzAnalysis::Init() {

  static ClassDocumentation<TopDalitzAnalysis> documentation
    ("There is no documentation for the TopDalitzAnalysis class");

}


tPVector TopDalitzAnalysis::particleID(PPtr top,tPVector final)
{
  PPtr orig=top;
  ShowerParticlePtr sorig;
  //////////////////////////////////////////////
  // Find the original top before it showers. //
  //////////////////////////////////////////////
  do
  {
      if(sqrt(orig->parents()[0]->momentum().m2())>=175001.) break;
      orig=orig->parents()[0];
      sorig=dynamic_ptr_cast<ShowerParticlePtr>(orig);
      if(sorig)
      {
	  if(sorig->perturbative()==2) break;
      }
  }
  while(dynamic_ptr_cast<ShowerParticlePtr>(orig));
  if(abs(orig->parents()[0]->id())!=ParticleID::t) orig=PPtr();
  // *top is always the thing that decays to bW. It's momentum is always 
  // exactly equal to that of the top when it is on-shell entering the 
  // decay shower, irrespective of whether there was (t)ISR or not! 
  // *orig:
  // If there is (t)ISR then *orig is the on-shell top right before it 
  // starts emitting gluons. In this case the momenta and masses of *top
  // and *orig are identical (e.g. evts 11,104,195,204,219).
  // *** WARNING(?) ***
  // If there is no (t)ISR then *orig is the top at production right
  // where it's mass is initially increased i.e. just prior to it's FSR
  // emission of gluons (i.e. where it's mass is maximally off-shell 
  // > 175 GeV). In this case the momenta and masses of *top and *orig
  // are different. If required this could simply be prevented by adding 
  // if(sqrt(orig->parents()[0]->momentum().m2())>=175001.) break;
  // ...or maybe even... 
  // if(sqrt(orig->parents()[0]->momentum().m2())>=175100.) break;
  // at the start of the above do-while loop (note: you should
  // check this mass cut is compatible with the top width i.e. this
  // assumes the top width is zero).


  ///////////////////////////////////////////////
  // Find the b and W before they shower/decay //
  ///////////////////////////////////////////////
  // Wboson turns out to be the W boson before the W boson 
  // which decays to ffbar. 
  // borig is the b quark which is one after the b quark which is 
  // produced at the point where the top decays. The b quark that comes
  // after borig splits to b+gluon.
  // The b quark and W boson that borig and Wboson
  // point to have exactly the same momentum as the b quark 
  // and W boson that come out of the top decay.
  // In the event that the hard matrix element correction is applied
  // then borig is equal to the momentum of the b after it has emitted 
  // the hard gluon.
  // This is a more straightforward bit of code than for the top (above).
  PPtr Wboson; 
  PPtr borig ;
  for(unsigned int ix=0;ix<top->children().size();++ix) {
      if(abs(top->children()[ix]->id())==ParticleID::Wplus) {
	  Wboson=top->children()[ix];
	  while(abs(Wboson->children()[0]->id())==ParticleID::Wplus&&
		!dynamic_ptr_cast<ShowerParticlePtr>(Wboson))
	      Wboson=Wboson->children()[0];
      } else {
	  borig=top->children()[ ix];
	  if(abs(borig->children()[0]->id())==ParticleID::b)
	      borig=borig->children()[0];
      }
  }
  // This bit fixes borig to point to the correct b-quark in the event
  // that a hard matrix element correction occurs.
  // If a hard MEC occurs you see a top decay to a bW then that b goes
  // *right away* to b+g instead of the normal b->b->b->b+g type chain.
  // Also if a hard MEC occurs the right W (for momentum conservation) 
  // is the W that decays, so we need to move Wboson down the family 
  // tree by one.
  if(borig->parents()[0]->children().size()>1) {
      if(abs(borig->parents()[0]->id())==ParticleID::b) 
	  { borig=borig->parents()[0]; }
  } 
  //////////////////////////////////////
  // Associate final states to t/W/b. //
  //////////////////////////////////////
  // The next bit tries to get borig to be the b-quark at the point where
  // it's scale is reset to begin FSR showering. This involves moving down
  // the tree by one.
  else if(borig->children()[0]->children().size()>1) {
      if((abs(borig->children()[0]->children()[0]->id())==ParticleID::b&&
	  abs(borig->children()[0]->children()[1]->id())==ParticleID::g)||
	 (abs(borig->children()[0]->children()[0]->id())==ParticleID::b&&
	  abs(borig->children()[0]->children()[1]->id())==ParticleID::g))
      {
          if(borig->children().size()==1)
	  borig=borig->children()[0];
      }
  }
  // The next for loop finds where the final state particles originated 
  // from i.e. b quarks or top quarks.
  tPVector tprod ;
  tPVector bottom;
  tPVector Wprod ;
  for(unsigned int ix=0;ix<final.size();++ix) {
      tPPtr part=final[ix]; 
      tPPtr last=part->parents()[0];
      tShowerParticlePtr shower;
      // This do-while hunts the parents of final state particle final[ix].
      do {
	  part=last;
          // If the particle has no parents set pointers to null & hang up.
	  if(part->parents().empty()) {last  =tPPtr();}
//	                               shower=tShowerParticlePtr();}
          // If we make it back to orig (top) break out.
          else if(last&&last==orig)   {break;}
          // If we make it back to Wboson (W boson) break out.
	  else if(last&&last==Wboson) {break;}
          // If we make it back to borig (bottom) break out.
	  else if(last&&last==borig)  {break;}
          // If none of the above set last to be the particle's parent,
          // and shower to be:
	  else {
	      last=last->parents()[0];
	      shower=dynamic_ptr_cast<ShowerParticlePtr>(last);
	  }
      } 
//      while(!part->parents().empty()&&shower);
      while(!part->parents().empty());
      // If the particle is found to originate from Wboson,borig,orig
      // then add it to the appropriate list.
      if(last) {
	  if(last==orig)        { tprod.push_back(final[ix]);  }
	  else if(last==borig)  { bottom.push_back(final[ix]); }
	  else if(last==Wboson) { Wprod.push_back(final[ix]);  } 
      }
  }
  // All of the final state particles (at the step 2 level) belonging
  // to the top/bottom/W boson have now been collected in tprod/
  // bottom/Wprod respectively.
  if(bottom.empty()) bottom.push_back(borig);
  if(Wprod.empty())  Wprod.push_back(Wboson);
  ////////////////////////
  // Check the momenta? //
  ////////////////////////
//  cout << "\n\ntprod.size()  " << tprod.size()  << endl;
//  cout << "bottom.size() " << bottom.size() << endl;
//   Lorentz5Momentum tq,bq,wb,diff;
//   for(unsigned int ix=0;ix<tprod.size() ;++ix) tq+=tprod[ix]->momentum();
//   for(unsigned int ix=0;ix<bottom.size();++ix) bq+=bottom[ix]->momentum();
//   for(unsigned int ix=0;ix<Wprod.size() ;++ix) wb+=Wprod[ix]->momentum();
//   cout << endl;
//   cout << tq << "\n" << bq << "\n" << wb << "\n"; 
//   diff = orig->momentum()-tq-bq-wb;
//   if(fabs(diff.t())>0.0001||fabs(diff.x())>0.0001||
//      fabs(diff.y())>0.0001||fabs(diff.z())>0.0001) 
//   { cout << "top quark   : " << diff  << endl; }
//   if(bq.m()>5000.1||bq.m()<4999.9) {
//       diff = bq-borig->momentum();
//   } else {
//       diff = bq-borig->children()[0]->momentum();
//   }
//   if(fabs(diff.t())>0.000000001||fabs(diff.x())>0.000000001||
//      fabs(diff.y())>0.000000001||fabs(diff.z())>0.000000001) 
//   { cout << "b jet       : " << diff   << endl; 
//     cout << "bq          : " << bq     << endl; 
//     cout << "*borig      : " << *borig << endl; 
//   }
//   diff = wb-Wboson->momentum();
//   if(fabs(diff.t())>0.00001||fabs(diff.x())>0.00001||
//      fabs(diff.y())>0.00001||fabs(diff.z())>0.00001) 
//   { cout << "W boson     : " << diff  << endl; }
//   cout << endl;
  ////////////////////
  // Jet Clustering //
  ////////////////////
  // First put all the products from the top and bottom quarks
  // into a single vector (bottom).
  for(unsigned int ix=0;ix<tprod.size();++ix)
    { bottom.push_back(tprod[ix]); }

  bottom.push_back(top); 
  bottom.push_back(orig); 
  return bottom;
}

void TopDalitzAnalysis::dalitz(tPVector finalPartons)
{
  PPtr top,orig;
  top  = finalPartons[finalPartons.size()-2];
  orig = finalPartons[finalPartons.size()-1];
  finalPartons.pop_back(); // pop out the on-shell top
  finalPartons.pop_back(); // pop out the top before it decays to bW
  ////////////////////
  // Jet Clustering //
  ////////////////////
  // For the Dalitz plot we need to have at least one gluon produced in 
  // the decay! Therefore we only try and plot a point if finalPartons.size()>1 
  // since the finalPartons array will only have one entry if there are no 
  // gluons radiated viz the b-quark which did not emit any other radiation
  // from the top or bottom makes finalPartons.size()>1 .
  if(finalPartons.size()>1)
    {
      _kint.clearMap();
      // Get KtJet to find two jets out of the list bottom.
      // Note, if the top did not radiate (tprod.size()==0) then KtJet
      // is giving back two jets made from the b and whatever the b
      // radiated. In this case the two jet momenta add to give bq. If
      // the top radiates then the two jet momenta seem to (and do) equal 
      // tq and bq respectively!
      KtJet::KtEvent ev = 
	  KtJet::KtEvent(_kint.convertToKtVectorList(finalPartons),1,1,1);
      ev.findJetsN(2);
      // Get the two jets ordered by their Pt (largest Pt first?).
      vector<KtJet::KtLorentzVector> ktjets = ev.getJetsPt();
      // Identify the jets.
      int nquark[2] = {0,0},iq;
      for(int ix = 0;ix<ev.getNConstituents();++ix)
	{
	  iq = ktjets[1].contains(*ev.getConstituents()[ix]);
          // iq = 1 if the 2nd jet (ktjets[1]) DOES contain constituent ix.
          // iq = 0 if the 2nd jet (ktjets[1]) DOESN'T contain constituent ix so
          // iq = 0 if the 1st jet (ktjets[0]) DOES contain constituent ix.
	  if(finalPartons[ix]->id() == -5) --nquark[iq];
          // If constituent ix is a bbar lower nquark[iq] by 1.
	  else if(finalPartons[ix]->id() == 5) ++nquark[iq];
          // If constituent ix is a b increase nquark[iq] by 1.
	}
      // Therefore nquark[0] is the number of b's-bbar's in jet 0  
      // Therefore nquark[1] is the number of b's-bbar's in jet 1. 
      Lorentz5Momentum pb,pg;
      // Identify the jets as being due to b/bbar quarks or gluons.
      if(top->id()>0)
	{
	  // If the decay is t->bW+ then if nquark[0] has more b-bbar than 
          // nquark[1] jet ktjet[0] must be the b jet.
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
	  // If the decay is tbar->bbarW- then if nquark[1] has less b-bbar than 
          // nquark[0] jet ktjet[1] must be the bbar jet.
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
      // Boost to the rest frame
      Hep3Vector boost;
      if(orig) boost=-orig->momentum().boostVector();
      else boost=-top->momentum().boostVector();
      pg.boost(boost);
      pb.boost(boost);
      Energy mt(top->mass());
      double xg(2.*pg.e()/mt),xb(2.*pb.e()/mt);
      if(_nout<50000) {
	_output[0] << xg << " " << 2.-xb-xg << "\n";
      }
      ++_nout;
    }
  return;
}

void TopDalitzAnalysis::threeJetAnalysis(Energy2 s,tPVector top, tPVector antitop)
{
  ////////////////////
  // Jet Clustering //
  ////////////////////
  // Chuck out the two top quarks *top and *orig first.
  top.pop_back();     // pop out the on-shell top 
  top.pop_back();     // pop out the top before it decays to bW
  antitop.pop_back(); 
  antitop.pop_back();
  tPVector finalPartons(top);
  // Now put the two lists of final state QCD particles together.
  for (unsigned int ix=0;ix<antitop.size();++ix) 
      finalPartons.push_back(antitop[ix]);
  // Cluster everything into 3 jets (if possible).
  if(finalPartons.size()>2)
    {
      _kint.clearMap();
      // Get KtJet to find two jets out of the list bottom.
      // Note, if the top did not radiate (tprod.size()==0) then KtJet
      // is giving back two jets made from the b and whatever the b
      // radiated. In this case the two jet momenta add to give bq. If
      // the top radiates then the two jet momenta seem to (and do) equal 
      // tq and bq respectively!
      KtJet::KtEvent ev = 
	  KtJet::KtEvent(_kint.convertToKtVectorList(finalPartons),1,1,1);
      ev.findJetsN(3);
      // Get the two jets ordered by their Pt (largest Pt first?).
      vector<KtJet::KtLorentzVector> ktjets = ev.getJetsPt();
//       // Identify the jets.
//       int nbs[3] = {0,0,0};
//       for(int ix = 0;ix<ev.getNConstituents();++ix)
//       {
// 	  for(unsigned int jx=0;jx<3;++jx) {
// 	      if(ktjets[jx].contains(*ev.getConstituents()[ix]))
// 	      { 
// 		  if(finalPartons[ix]->id() == -5) --nbs[jx];
// 		  if(finalPartons[ix]->id() ==  5) ++nbs[jx];
// 	      }
// 	  }
//       }
      // Therefore nbs[ix] is the number of b's-bbar's in jet ix.  
      Lorentz5Momentum p0(ktjets[0]),p1(ktjets[1]),p2(ktjets[2]);

      ///////////////////////
      // Calculate delta R //
      ///////////////////////
      double deltaR(0.);
      deltaR =            sqrt(sqr(p0.eta()-p1.eta())+sqr(p0.phi()-p1.phi())) ;
      deltaR = min(deltaR,sqrt(sqr(p0.eta()-p2.eta())+sqr(p0.phi()-p2.phi())));
      deltaR = min(deltaR,sqrt(sqr(p1.eta()-p2.eta())+sqr(p1.phi()-p2.phi())));
      // If jets pass Et and deltaR separation cuts then add a 
      // point to the histogram.
      if(deltaR>0.7&&p0.et()>10000.&&p1.et()>10000.&&p2.et()>10000.) _deltaR += deltaR;

      //////////////////
      // Calculate y3 //
      //////////////////
      double y3(0.);
      Hep3Vector np0,np1,np2;
      np0    = p0.vect()/p0.vect().mag(); 
      np1    = p1.vect()/p1.vect().mag(); 
      np2    = p2.vect()/p2.vect().mag(); 
      y3 =        (2./s)*min(sqr(p0.e()),sqr(p1.e()))*(1.-np0*np1) ;
      y3 = min(y3,(2./s)*min(sqr(p0.e()),sqr(p2.e()))*(1.-np0*np2));
      y3 = min(y3,(2./s)*min(sqr(p1.e()),sqr(p2.e()))*(1.-np1*np2));
      // If jets pass Et and deltaR separation cuts then add a 
      // point to the histogram.
/// Do log to base 10 ?
// Yes:
      if(deltaR>0.7&&p0.et()>10000.&&p1.et()>10000.&&p2.et()>10000.) _logy3 += log(y3)/log(10.);
// No:
//      if(deltaR>0.7&&p0.et()>10000.&&p1.et()>10000.&&p2>10000.) _logy3 += log(y3);
    }
  return;
}

void TopDalitzAnalysis::dofinish() {
  AnalysisHandler::dofinish();

  /////////////////
  // Dalitz Plot //
  /////////////////
  _output[0] << "PLOT\n";
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
    {_output[0] << xgg[ix] << " " << upper[ix] << "\n";}
  for(int ix=lower.size()-1;ix>=0;--ix)
    {_output[0] << xgg[ix] << " " << lower[ix] << "\n";}
  _output[0] << "JOIN RED " << "\n";
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
      if(xg<xgmax) _output[0] << xg << " " << xa << "\n";
    }
  _output[0] << "JOIN BLUE" << endl;
  // phase space for radiation from top
  kappa=1+0.25*sqr(1.-a+c+lam)/(kappa-c);
  cerr << "top    kappa " << kappa << endl;
  double u,w,v;
  //double zmin=1.-(1.-a)/(kappa+2.*sqrt(a*(kappa-1.)))+0.00001;
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
      if(xg<xgmax) _output[0] << xg << " " << xa << "\n";
    }
  _output[0] << "JOIN GREEN" << endl;
  _output[0].close();

  ////////////
  // DeltaR //
  ////////////
  _deltaR.topdrawOutput(_output[1],true,false,false,false,"RED","delta(R)");
  _output[1].close();

  /////////////
  // log(y3) //
  /////////////
  _logy3.topdrawOutput(_output[2],true,false,false,false,"RED","log(y3)");
  _output[2].close();

}
