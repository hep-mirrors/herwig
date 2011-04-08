// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2CharginoNeutralinoPowheg class.
//

#include "MEPP2CharginoNeutralinoPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/Susy/SusyBase.h"
#include "Herwig++/MatrixElement/HardVertex.h"
#include <numeric>

using namespace Herwig;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;


MEPP2CharginoNeutralinoPowheg::MEPP2CharginoNeutralinoPowheg() 
  : process_(0), maxFlavour_(5) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

void MEPP2CharginoNeutralinoPowheg::doinit() {
  MEPP2GauginoGauginoPowheg::doinit();
  // get the photon and W ParticleData objects
  Wplus_    = getParticleData(ThePEG::ParticleID::Wplus);
  Wminus_   = getParticleData(ThePEG::ParticleID::Wminus);

  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2CharginoNeutralinoPowheg::doinit" 
			  << Exception::abortnow;
  // do the initialisation (see Herwig::SusyBase Class)
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  CNWVertex_ = hwsm->vertexCNW();
  CFSVertex_ = hwsm->vertexCFSF();
  NFSVertex_ = hwsm->vertexNFSF();
  GSSVertex_ = hwsm->vertexGSFSF();

}

Selector<const ColourLines *>
MEPP2CharginoNeutralinoPowheg::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2"), c2("1 2 -3");
  Selector<const ColourLines *> sel;
  if(abs(diag->id())==1)
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEPP2CharginoNeutralinoPowheg::getDiagrams() const {
  // loop over the processes we need
  tcPDPtr cha[2] = {getParticleData(1000024),getParticleData(1000037)};
  tcPDPtr chab[2];
  tcPDPtr neu[4] = {getParticleData(1000022),getParticleData(1000023),
		    getParticleData(1000025),getParticleData(1000035)};
  tcPDPtr q,qb;
  for(unsigned int x=0;x<2;++x)
    chab[x] = cha[x]->CC();

  for(int i = 1; i <= maxFlavour_; ++i) {
    if(i%2==1){
      q  = getParticleData(i);
      qb = getParticleData(i+1)->CC();}
    else{
      q  = getParticleData(i);
      qb = getParticleData(i-1)->CC();}

    tcPDPtr qLt = getParticleData(1000000+i);
    tcPDPtr qRt = getParticleData(2000000+i);
    tcPDPtr qLu = getParticleData(1000000+i+1);
    tcPDPtr qRu = getParticleData(2000000+i+1);
    tcPDPtr qLd = getParticleData(1000000+i-1);
    tcPDPtr qRd = getParticleData(2000000+i-1);
    for(unsigned int ix=0;ix<2;++ix){
      for(unsigned int jx=0;jx<4;++jx){
	// Wminus mediated processes
	if(process_==-1 || process_==4*ix+jx+1){
	  // q is down type
	  if(i%2==1) {
	    // W-mediated s-channel
	    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wminus_,
			 3, chab[ix], 3, neu[jx], -1)));
	    
	    // ~qL mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qLt, qb,
			 3, chab[ix], 1, neu[jx], -2)));
	    // ~qR mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qRt, qb,
			 3, chab[ix], 1, neu[jx], -3)));
	    
	    // ~qL mediated u-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qLu, qb,
			 3, neu[jx], 1, chab[ix], -4)));
	    // ~qR mediated u-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qRu, qb,
			 3, neu[jx], 1, chab[ix], -5)));
	  }
	}
	// Wplus mediated processes
	if(process_ == 0 || process_ == 4*ix+jx+9){
	  // q is up type
	  if(i%2==0) {
	    // W-mediated s-channel
	    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Wplus_,
			 3, cha[ix], 3, neu[jx], -1)));

	    // ~qL mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qLt, qb,
			 3, cha[ix], 1, neu[jx], -2)));
	    // ~qR mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qRt, qb,
			 3, cha[ix], 1, neu[jx], -3)));

	    // ~qL mediated u-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qLu, qb,
			 3, neu[jx], 1, cha[ix], -4)));
	    // ~qR mediated u-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qRu, qb,
			 3, neu[jx], 1, cha[ix], -5)));
	  }
	}
      }
    }
  }
}

unsigned int MEPP2CharginoNeutralinoPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2CharginoNeutralinoPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2CharginoNeutralinoPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2CharginoNeutralinoPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2CharginoNeutralinoPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFWVertex_ << FFGVertex_ << CNWVertex_ << CFSVertex_ << NFSVertex_
     << GSSVertex_ << Wplus_ << Wminus_ << process_ << maxFlavour_;
}

void MEPP2CharginoNeutralinoPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> FFGVertex_ >> CNWVertex_ >> CFSVertex_ >> NFSVertex_
     >> GSSVertex_ >> Wplus_ >> Wminus_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2CharginoNeutralinoPowheg> 
MEPP2CharginoNeutralinoPowheg::initMEPP2CharginoNeutralinoPowheg;
// Definition of the static class description member.

void MEPP2CharginoNeutralinoPowheg::Init() {

  static ClassDocumentation<MEPP2CharginoNeutralinoPowheg> documentation
    ("MEPP2CharginoNeutralinoPowheg implements the ME calculation"
     " of the fermion-antifermion to chargino-neutralino"
     " hard process.");

  static Switch<MEPP2CharginoNeutralinoPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2CharginoNeutralinoPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAllWm
    (interfaceProcess,
     "All",
     "Generate all Wminus-mediated processes"
     " (i.e. all combinations of chargino"
     "and neutralino mass eigenstate pairs.)",
     -1);
  static SwitchOption interfaceProcessAllWp
    (interfaceProcess,
     "All",
     "Generate all Wplus-mediated processes"
     " (i.e. all combinations of chargino"
     "and neutralino mass eigenstate pairs.)",
     0);
  static SwitchOption interfaceProcessCharginom11
    (interfaceProcess,
     "chi-1chi01",
     "Only produce chi-1, chi01 pairs.",
     1);
  static SwitchOption interfaceProcessCharginom12
    (interfaceProcess,
     "chi-1chi02",
     "Only produce chi-1, chi02 pairs.",
     2);
  static SwitchOption interfaceProcessCharginom13
    (interfaceProcess,
     "chi-1chi03",
     "Only produce chi-1, chi03 pairs.",
     3);
  static SwitchOption interfaceProcessCharginom14
    (interfaceProcess,
     "chi-1chi04",
     "Only produce chi-1, chi04 pairs.",
     4);
  static SwitchOption interfaceProcessCharginom21
    (interfaceProcess,
     "chi-2chi01",
     "Only produce chi-2, chi01 pairs.",
     5);
  static SwitchOption interfaceProcessCharginom22
    (interfaceProcess,
     "chi-2chi02",
     "Only produce chi-2, chi02 pairs.",
     6);
  static SwitchOption interfaceProcessCharginom23
    (interfaceProcess,
     "chi-2chi03",
     "Only produce chi-2, chi03 pairs.",
     7);
  static SwitchOption interfaceProcessCharginom24
    (interfaceProcess,
     "chi-2chi04",
     "Only produce chi-2, chi04 pairs.",
     8);
  static SwitchOption interfaceProcessCharginop11
    (interfaceProcess,
     "chi+1chi01",
     "Only produce chi+1, chi01 pairs.",
     9);
  static SwitchOption interfaceProcessCharginop12
    (interfaceProcess,
     "chi+1chi02",
     "Only produce chi+1, chi02 pairs.",
     10);
  static SwitchOption interfaceProcessCharginop13
    (interfaceProcess,
     "chi+1chi03",
     "Only produce chi+1, chi03 pairs.",
     11);
  static SwitchOption interfaceProcessCharginop14
    (interfaceProcess,
     "chi+1chi04",
     "Only produce chi+1, chi04 pairs.",
     12);
  static SwitchOption interfaceProcessCharginop21
    (interfaceProcess,
     "chi+2chi01",
     "Only produce chi+2, chi01 pairs.",
     13);
  static SwitchOption interfaceProcessCharginop22
    (interfaceProcess,
     "chi+2chi02",
     "Only produce chi+2, chi02 pairs.",
     14);
  static SwitchOption interfaceProcessCharginop23
    (interfaceProcess,
     "chi+2chi03",
     "Only produce chi+2, chi03 pairs.",
     15);
  static SwitchOption interfaceProcessCharginop24
    (interfaceProcess,
     "chi+2chi04",
     "Only produce chi+2, chi04 pairs.",
     16);


  static Parameter<MEPP2CharginoNeutralinoPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2CharginoNeutralinoPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

}

Selector<MEBase::DiagramIndex>
MEPP2CharginoNeutralinoPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1)     sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[2], i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2CharginoNeutralinoPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  output.finite =-8.+sqr(Constants::pi);
  output.finite *= loWeight();
  return output;
}

double MEPP2CharginoNeutralinoPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	vector<SpinorWaveFunction>    & spout ,
	vector<SpinorBarWaveFunction> & sbarout ,
	//ScalarWaveFunction & interqL,ScalarWaveFunction & interqR,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  tcPDPtr squark[4];
  if (abs(mePartonData()[0]->id())%2==0) {
    // u-channel
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())-1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())-1);
    // t-channel
    squark[2] = getParticleData(1000000+abs(mePartonData()[0]->id())  );
    squark[3] = getParticleData(2000000+abs(mePartonData()[0]->id())  );
  }
  else {
    // u-channel
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())+1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())+1);
    // t-channel
    squark[2] = getParticleData(1000000+abs(mePartonData()[0]->id())  );
    squark[3] = getParticleData(2000000+abs(mePartonData()[0]->id())  );
  }
  // conjugate spinors for t-channel diagram
  vector<SpinorWaveFunction> sbaroutconj;
  vector<SpinorBarWaveFunction> spoutconj;
  for(unsigned int of1=0;of1<2;++of1) {
    sbaroutconj.push_back(SpinorWaveFunction (-sbarout[of1].momentum(),
					      sbarout[of1].particle(),
					      sbarout[of1].wave().bar().conjugate(),
				 	      sbarout[of1].direction()));
    spoutconj.push_back(SpinorBarWaveFunction(-spout[of1].momentum(),
					      spout[of1].particle(),
					      spout[of1].wave().bar().conjugate(),
					      spout[of1].direction()));
  }
  // storage of the matrix elements for specific diagrams
  vector<double> me(5, 0.);
  double me2(0.);
  // storage of the individual diagrams
  vector<Complex> diag(5, 0.);
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin1Half, PDT::Spin1Half);
//   unsigned int propopt;
//   if(status()==RealQG || status()==RealQBarG){
//     propopt = 7;
//   }
//   else propopt = 3;
  unsigned int propopt = 3;
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interW = FFWVertex_->
	evaluate(q2, 1, mePartonData()[0]->iCharge() > 0 ? Wplus_ : Wminus_,
		 sp[if1],sbar[if2]);
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  // s-channel
	  diag[0] = CNWVertex_->evaluate(q2, spout[of1],  sbarout[of2], interW);
	  // t and u channels
	  for(unsigned int iq=0;iq<2;++iq) {
	    // t-channel
	    ScalarWaveFunction intersq = CFSVertex_->
	      evaluate(q2, propopt, squark [iq   ], sp[if1], spoutconj[of1]);
	    diag[iq+1] = -NFSVertex_->
	      evaluate(q2, sbaroutconj[of2], sbar[if2], intersq);
	    // u-channel
	    intersq = NFSVertex_->
		evaluate(q2, propopt, squark[iq+2], sp[if1], sbarout[of2]);
	    diag[iq+3] = CFSVertex_->
	      evaluate(q2, spout[of1], sbar[if2], intersq);
	  }
	  // Individual diagrams
	  for(unsigned int id=0;id<5;++id) me[id] += norm(diag[id]);
	  // sum up the matrix elements
	  Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  me2 += norm(total);
	  pme(if1, if2, of1, of2) = total;
	}
      }
    }
  }
  if(first) {
    DVector save(3);
    for(DVector::size_type ix = 0; ix < 3; ++ix)
      save[ix] = me[ix]/12.;
    meInfo(save);
    me_.reset(pme);
  }
  return me2/12.;
}

double MEPP2CharginoNeutralinoPowheg::loME(const cPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta,
				    bool first) const {
  // wavefunctions for the incoming fermions
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    sp[i]   = SpinorWaveFunction   (momenta[0], particles[0], i,
				    incoming);
    sbar[i] = SpinorBarWaveFunction(momenta[1], particles[1], i,
				    incoming);
  }
  // outgoing neutralino wavefunctions
  vector<SpinorWaveFunction> spout(2);
  vector<SpinorBarWaveFunction> sbarout(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    spout[i] = SpinorWaveFunction(momenta[2], particles[2], i,
				  outgoing);
    sbarout[i] = SpinorBarWaveFunction(momenta[3], particles[3], i,
				       outgoing);
  }
  return qqbarME(sp,sbar,spout,sbarout,first);
}


double MEPP2CharginoNeutralinoPowheg::realME(const cPDVector & particles,
 				      const vector<Lorentz5Momentum> & momenta) const {
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // squarks
  tcPDPtr squark[4];
  if (abs(mePartonData()[0]->id())%2==0) {
    // u-channel
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())-1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())-1);
    // t-channel
    squark[2] = getParticleData(1000000+abs(mePartonData()[0]->id())  );
    squark[3] = getParticleData(2000000+abs(mePartonData()[0]->id())  );
  }
  else {
    // u-channel
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())+1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())+1);
    // t-channel
    squark[2] = getParticleData(1000000+abs(mePartonData()[0]->id())  );
    squark[3] = getParticleData(2000000+abs(mePartonData()[0]->id())  );
  }

  // wavefunctions for the q qbar -> chi chi g process
  // q qbar
  if(abs(particles[0]->id())<=6 &&
     abs(particles[1]->id())<=6) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
    }
  }
  else if(particles[0]->id()>0 &&
	  particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      gluon[i]= VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  // g qbar
  else if(particles[0]->id()==ParticleID::g &&
	  particles[1]->id()<0) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      cerr << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the outgoing charginos/neutralinos
  vector<SpinorWaveFunction> spout(2),sbaroutconj(2);
  vector<SpinorBarWaveFunction> sbarout(2),spoutconj(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    spout[i]       = SpinorWaveFunction(momenta[2], particles[2], i,
					outgoing);
    spoutconj[i]   = SpinorBarWaveFunction(momenta[2], particles[2],
					   spout[i].wave().bar().conjugate(),
					   outgoing);
    sbarout[i]     = SpinorBarWaveFunction(momenta[3], particles[3], i,
					   outgoing);
    sbaroutconj[i] = SpinorWaveFunction   (momenta[3], particles[3],
					   sbarout[i].wave().bar().conjugate(),
					   outgoing);
  }
  double output(0.);
  vector<Complex> diag(14,0.);
  Energy2 q2 = scale();
  unsigned int propopt = 3;
//   if(status()==RealQG || status()==RealQBarG){
//     propopt = 7;
//   }
//   else propopt = 3;
  ScalarWaveFunction intersq,intersq2;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	SpinorWaveFunction    inters = FFGVertex_->evaluate(q2, 5, sp[ihel1].particle()->CC(),
							    sp[ihel1], gluon[ohel1]);
	SpinorBarWaveFunction interb = FFGVertex_->evaluate(q2,5,sbar[ihel2].particle()->CC(),
							    sbar[ihel2],gluon[ohel1]);	
	VectorWaveFunction interW[2] = {FFWVertex_->
					evaluate(q2, 1, mePartonData()[0]->iCharge() > 0 ?
						 Wplus_ : Wminus_, inters, sbar[ihel2]),
					FFWVertex_->
					evaluate(q2, 1, mePartonData()[0]->iCharge() > 0 ?
						 Wplus_ : Wminus_, sp[ihel1], interb)};
	
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  for(unsigned int ohel3=0;ohel3<2;++ohel3) {
	    // s-channel diagrams
	    // first W diagram (emission from quark)
	    diag[0] = CNWVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interW[0]);
	    // second W diagram (emission from anti-quark)
	    diag[1] = CNWVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interW[1]);
	    // t and u channels
	  for(unsigned int iq=0;iq<2;++iq) {
	    // t-channel
	    // emission from quark
	    intersq = CFSVertex_->
	      evaluate(q2, propopt, squark [iq], inters, spoutconj[ohel2]);
	    diag[6*iq+2] = -NFSVertex_->
	      evaluate(q2, sbaroutconj[ohel3], sbar[ihel2], intersq);
	    // emission from anti-quark
	    intersq = CFSVertex_->
	      evaluate(q2, propopt, squark [iq], sp[ihel1], spoutconj[ohel2]);
	    diag[6*iq+3] = -NFSVertex_->
	      evaluate(q2, sbaroutconj[ohel3], interb, intersq);
	    // emission from intermediate
	    intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq],gluon[ohel1],intersq);
	    diag[6*iq+4] = -NFSVertex_->
	      evaluate(q2,sbaroutconj[ohel3], sbar[ihel2], intersq2);

	    // u-channel
	    // emission from quark
	    intersq = NFSVertex_->
	      evaluate(q2, propopt, squark [iq+2], inters, sbarout[ohel3]);
	    diag[6*iq+5] = CFSVertex_->
	      evaluate(q2, spout[ohel2], sbar[ihel2], intersq);
	    // emission from anti-quark
	    intersq = NFSVertex_->
	      evaluate(q2, propopt, squark [iq+2], sp[ihel1], sbarout[ohel3]);
	    diag[6*iq+6] = CFSVertex_->
	      evaluate(q2, spout[ohel2], interb, intersq);
	    // emission from intermediate
	    intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq+2],gluon[ohel1],intersq);
	    diag[6*iq+7] = CFSVertex_->
	      evaluate(q2,spout[ohel2], sbar[ihel2], intersq2);
	  }


	    // add them up
	    Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    output += norm(total);

	    // Gauge invariance check
// 	    double all = 0.;
// 	    for(unsigned int ix=0;ix<14;++ix){
// 	      all += norm(diag[ix]);}
// 	    if(isnan(norm(total)/all) != 1)
// 	      cout << norm(total)/all << endl;


	  }
	}
      }
    }
  }
  // colour and spin factors
  if(particles[0]->id()==ParticleID::g ||
     particles[1]->id()==ParticleID::g){
    output *= 1./24.;
  }
  else {
    output *=1./9.;
  }
  // divided by 2 g_S^2
  return 0.5*output/norm(FFGVertex_->norm());
}

void MEPP2CharginoNeutralinoPowheg::constructVertex(tSubProPtr sub) {
//   //get particles
//   ParticleVector ext(4);
//   ext[0] = sub->incoming().first;
//   ext[1] = sub->incoming().second;
//   ext[2] = sub->outgoing()[0];
//   ext[3] = sub->outgoing()[1];
//   if( ext[0]->id() != mePartonData()[0]->id() ) swap(ext[0], ext[1]);
//   if( ext[2]->id() != mePartonData()[2]->id() ) swap(ext[2], ext[3]);

//   //First calculate wave functions with off-shell momenta
//   //to calculate correct spin information
//   vector<SpinorWaveFunction> sp;
//   SpinorWaveFunction(sp, ext[0], incoming, false);
//   vector<SpinorBarWaveFunction> sbar;
//   SpinorBarWaveFunction(sbar, ext[1], incoming, false);
//   vector<SpinorWaveFunction> spout;
//   SpinorWaveFunction(spout, ext[2], outgoing, false);
//   vector<SpinorBarWaveFunction> sbarout;
//   SpinorBarWaveFunction(sbarout, ext[3], outgoing, false);

//   //Need to use rescale momenta to calculate matrix element
//   cPDVector data(4);
//   vector<Lorentz5Momentum> momenta(4);
//   for( size_t i = 0; i < 4; ++i ) {
//     data[i] = ext[i]->dataPtr();
//     momenta[i] = ext[i]->momentum();
//   }
//   rescaleMomenta(momenta, data);
//   SpinorWaveFunction spr(rescaledMomenta()[0], data[0], incoming);
//   SpinorBarWaveFunction sbr(rescaledMomenta()[1], data[1],incoming);
//   SpinorWaveFunction spout1(rescaledMomenta()[2], data[2], outgoing);
//   SpinorBarWaveFunction sbarout2(rescaledMomenta()[3], data[3], outgoing);
//   for( unsigned int ihel = 0; ihel < 2; ++ihel ) {  
//     spr.reset(ihel);
//     sp[ihel] = spr;
//     sbr.reset(ihel);
//     sbar[ihel] = sbr;
//     spout1.reset(ihel);
//     spout[ihel] = spout1;
//     sbarout2.reset(ihel);
//     sbarout[ihel] = sbarout2;
//   }


//   qqbarME(sp, sbar, spout, sbarout, true);
//   HardVertexPtr hv = new_ptr(HardVertex());
//   hv->ME(me_);
//   for(unsigned int i = 0; i < 4; ++i )
//     dynamic_ptr_cast<SpinfoPtr>(ext[i]->spinInfo())->setProductionVertex(hv);  
}





