// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2CharginoCharginoPowheg class.
//

#include "MEPP2CharginoCharginoPowheg.h"
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


MEPP2CharginoCharginoPowheg::MEPP2CharginoCharginoPowheg() 
  : process_(0), maxFlavour_(5) {
  vector<unsigned int> mopt(2,1);
  massOption(mopt);
}

void MEPP2CharginoCharginoPowheg::doinit() {
  NLODrellYanBase::doinit();
  // get the photon and Z ParticleData objects
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2CharginoCharginoPowheg::doinit" 
			  << Exception::abortnow;
  // do the initialisation (see Herwig::SusyBase Class)
  FFZVertex_ = hwsm->vertexFFZ();
  FFGVertex_ = hwsm->vertexFFG();
  NNZVertex_ = hwsm->vertexNNZ();
  NFSVertex_ = hwsm->vertexNFSF();
}

Selector<const ColourLines *>
MEPP2CharginoCharginoPowheg::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2"), c2("1 2 -3");
  Selector<const ColourLines *> sel;
  if(abs(diag->id())==1)
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEPP2CharginoCharginoPowheg::getDiagrams() const {
  // loop over the processes we need
  tcPDPtr chi[4] = {getParticleData(1000022),getParticleData(1000023),
		    getParticleData(1000025),getParticleData(1000035)};
  for(int i = 1; i <= maxFlavour_; ++i) {
    tcPDPtr q  = getParticleData(i);
    tcPDPtr qb = q->CC();
    tcPDPtr qL = getParticleData(1000000+i);
    tcPDPtr qR = getParticleData(2000000+i);
    for(int c1=0;c1<4;++c1) {
      for(int c2=0;c2<=c1;++c2) {
	if(process_==0 || process_ == 4*c2+c1+1) {
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_   ,
		       3, chi[c1], 3, chi[c2], -1)));
	  add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
		       1, chi[c1], 3, chi[c2], -2)));
	  add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
		       1, chi[c1], 3, chi[c2], -3)));
	  add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
		       3, chi[c1], 1, chi[c2], -4)));
	  add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
		       3, chi[c1], 1, chi[c2], -5)));

	  //cout << c1 << " " << c2 << " " << 4*c2+c1+1 << endl;

	}
      }
    }
  }
}

unsigned int MEPP2CharginoCharginoPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2CharginoCharginoPowheg::orderInAlphaEW() const {
  return 2;
}

IBPtr MEPP2CharginoCharginoPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2CharginoCharginoPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2CharginoCharginoPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFGVertex_ << NNZVertex_ << NFSVertex_
     << Z0_ << process_ << maxFlavour_;
}

void MEPP2CharginoCharginoPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFGVertex_ >> NNZVertex_ >> NFSVertex_
     >> Z0_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2CharginoCharginoPowheg> 
MEPP2CharginoCharginoPowheg::initMEPP2CharginoCharginoPowheg;
// Definition of the static class description member.

void MEPP2CharginoCharginoPowheg::Init() {

  static ClassDocumentation<MEPP2CharginoCharginoPowheg> documentation
    ("MEPP2CharginoCharginoPowheg implements the ME calculation"
     " of the fermion-antifermion to chargino-chargino or"
     " neutralino-neutralino hard process.");


  static Switch<MEPP2CharginoCharginoPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2CharginoCharginoPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all the processes"
     " (i.e. both chargino and neutralino pair production"
     " of all mass eigenstates.)",
     0);
  static SwitchOption interfaceProcessNeutralino11
    (interfaceProcess,
     "chi01chi01",
     "Only produce chi01, chi01 pairs.",
     1);
  static SwitchOption interfaceProcessNeutralino12
    (interfaceProcess,
     "chi01chi02",
     "Only produce chi01, chi02 pairs.",
     2);
  static SwitchOption interfaceProcessNeutralino13
    (interfaceProcess,
     "chi01chi03",
     "Only produce chi01, chi03 pairs.",
     3);
  static SwitchOption interfaceProcessNeutralino14
    (interfaceProcess,
     "chi01chi04",
     "Only produce chi01, chi04 pairs.",
     4);
  static SwitchOption interfaceProcessNeutralino22
    (interfaceProcess,
     "chi02chi02",
     "Only produce chi02, chi02 pairs.",
     6);
  static SwitchOption interfaceProcessNeutralino23
    (interfaceProcess,
     "chi02chi03",
     "Only produce chi02, chi03 pairs.",
     7);
  static SwitchOption interfaceProcessNeutralino24
    (interfaceProcess,
     "chi02chi04",
     "Only produce chi02, chi04 pairs.",
     8);
  static SwitchOption interfaceProcessNeutralino33
    (interfaceProcess,
     "chi03chi03",
     "Only produce chi03, chi04 pairs.",
     11);
  static SwitchOption interfaceProcessNeutralino34
    (interfaceProcess,
     "chi03chi04",
     "Only produce chi03, chi04 pairs.",
     12);
  static SwitchOption interfaceProcessNeutralino44
    (interfaceProcess,
     "chi04chi04",
     "Only produce chi04, chi04 pairs.",
     16);

  static Parameter<MEPP2CharginoCharginoPowheg,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The maximum flavour of the incoming quarks",
     &MEPP2CharginoCharginoPowheg::maxFlavour_, 5, 1, 5,
     false, false, Interface::limited);

}


Selector<MEBase::DiagramIndex>
MEPP2CharginoCharginoPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 )      sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[2], i);
    else if ( diags[i]->id() == -4 ) sel.insert(meInfo()[3], i);
    else if ( diags[i]->id() == -5 ) sel.insert(meInfo()[4], i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2CharginoCharginoPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  output.finite =-8.+sqr(Constants::pi);
  return output;
}

// ofstream myfile ("Our_MSSM_scale.txt");
// if (myfile.is_open())
//   {
//     myfile << "Our scale " << scale() << endl;
//     myfile.close();
//   }
// else cout << "Unable to open file";

double MEPP2CharginoCharginoPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	vector<SpinorWaveFunction>    & spout ,
	vector<SpinorBarWaveFunction> & sbarout ,
	//ScalarWaveFunction & interqL,ScalarWaveFunction & interqR,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  // squarks for the t-channel
  tcPDPtr squark[2]= {getParticleData(1000000+abs(mePartonData()[0]->id())),
		      getParticleData(2000000+abs(mePartonData()[0]->id()))};
  // conjugate spinors for t-channel exchange diagram
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
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interV = FFZVertex_->
	evaluate(q2, 1, Z0_, sp[if1],sbar[if2]);
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  // s-channel Z exchange
	  diag[0] = NNZVertex_->evaluate(q2, spout[of1],  sbarout[of2], interV);
	  // t-channel squark exchanges	  
	  for(unsigned int iq=0;iq<2;++iq) {
	    // 1st t-channel
	    ScalarWaveFunction intersq = NFSVertex_->
	      evaluate(q2, 3, squark[iq], sp[if1], sbarout[of2]);
	    diag[2*iq+1] = 
	      NFSVertex_->evaluate(q2, spout[of1], sbar[if2], intersq);
	    // swapped t-channel
	    intersq = NFSVertex_->
	      evaluate(q2, 3, squark[iq], sp[if1], spoutconj[of1]);
	    diag[2*iq+2] = 
	      -NFSVertex_->evaluate(q2, sbaroutconj[of2], sbar[if2], intersq);
	  }
	  // individual diagrams
	  for(unsigned int id=0;id<5;++id){
	    me[id] += norm(diag[id]);
	  }
	  // sum up the matrix elements
	  Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  me2 += norm(total);
	  pme(if1, if2, of1, of2) = total;
	}
      }
    }
  }
  if(first) {
    DVector save(5);
    for(DVector::size_type ix = 0; ix < 5; ++ix)
      save[ix] = me[ix]/12.;
    meInfo(save);
    me_.reset(pme);
  }
  if(mePartonData()[2]->id() == mePartonData()[3]->id()) me2 *= 0.5;
  return me2/12.;
}

double MEPP2CharginoCharginoPowheg::loME(const cPDVector & particles,
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


double MEPP2CharginoCharginoPowheg::realME(const cPDVector & particles,
 				      const vector<Lorentz5Momentum> & momenta) const {
//   vector<SpinorWaveFunction> sp(2);
//   vector<SpinorBarWaveFunction> sbar(2);
//   vector<VectorWaveFunction> gluon(2);
//   // wavefunctions for the q qbar -> sf sf g process
//   if(particles[0]->id()==-particles[1]->id()) {
//     for( unsigned int i = 0; i < 2; ++i ) {
//       sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
//       sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
//       gluon[i]= VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
//     }
//   }
//   else if(particles[0]->id()==ParticleID::g &&
// 	  particles[1]->id()<0) {
//     for( unsigned int i = 0; i < 2; ++i ) {
//       sp[i]   = SpinorWaveFunction   (momenta[4],particles[4],  i,outgoing);
//       sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
//       gluon[i]= VectorWaveFunction   (momenta[0],particles[0],2*i,incoming);
//     }
//   }
//   else if(particles[0]->id()>0 &&
// 	  particles[1]->id()==ParticleID::g) {
//     for( unsigned int i = 0; i < 2; ++i ) {
//       sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
//       sbar[i] = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
//       gluon[i]= VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
//     }
//   }
//   else {
//     for(unsigned int ix=0;ix<particles.size();++ix) {
//       cerr << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
//     }
//     assert(false);
//   }
//   // wavefunctions for the outgoing charginos
//   vector<SpinorWaveFunction> spout(2);
//   vector<SpinorBarWaveFunction> sbarout(2);
//   for( unsigned int i = 0; i < 2; ++i ) {
//     spout[i] = SpinorWaveFunction(momenta[2], particles[2], i,
// 			       incoming);
//     sbarout[i] = SpinorBarWaveFunction(momenta[3], particles[3], i,
// 				    incoming);
//   }

//   double output(0.);
//   Complex diag[4]={0.,0.,0.,0.};
//   Energy2 shat = scale();
//   for(unsigned int ihel1=0;ihel1<2;++ihel1) {
//     for(unsigned int ihel2=0;ihel2<2;++ihel2) {
//       for(unsigned int ohel1=0;ohel1<2;++ohel1) {
// 	// first Z diagram (emission from quark)
//  	SpinorWaveFunction inters = FFGVertex_->evaluate(shat, 5, sp[ihel1].particle(),
// 							 sp[ihel1], gluon[ohel1]);
// 	VectorWaveFunction interV = FFZVertex_->evaluate(shat, 1, Z0_, inters, 
// 							 sbar[ihel2]);
// 	diag[0] = CCZVertex_->evaluate(shat, spout[ihel1], sbarout[ihel2], interV);
// 	// second Z diagram (emission from anti-quark)
// 	SpinorBarWaveFunction interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
// 							    sbar[ihel1],gluon[ohel1]);
// 	interV = FFZVertex_->evaluate(shat, 1, Z0_, sp[ihel1], 
// 				      interb);
// 	diag[1] = CCZVertex_->evaluate(shat, spout[ihel1], sbarout[ihel2], interV);


// 	if(particles[2]->id()==-particles[3]->id()&&particles[2]->charged()) {
//  	  // first photon diagram
//  	inters = FFGVertex_->evaluate(shat, 5, sp[ihel1].particle(),
// 							 sp[ihel1], gluon[ohel1]);
// 	interV = FFPVertex_->evaluate(shat, 1, gamma_, inters, 
// 							 sbar[ihel2]);
// 	diag[2] = CCZVertex_->evaluate(shat, spout[ihel1], sbarout[ihel2], interV);
// 	// second photon diagram
// 	interb = FFGVertex_->evaluate(shat,5,sbar[ihel1].particle(),
// 							    sbar[ihel1],gluon[ohel1]);
// 	interV = FFPVertex_->evaluate(shat, 1, gamma_, sp[ihel1], 
// 				      interb);
// 	diag[3] = CCZVertex_->evaluate(shat, spout[ihel1], sbarout[ihel2], interV);
//  	}
// 	// add them up
// 	output += norm(diag[0]+diag[1]+diag[2]+diag[3]);
//       }
//     }
//   }
//   // colour and spin factors
//   if(particles[0]->id()==-particles[1]->id()) {
//     output *= 1./9.;
//   }
//   else  {
//     output *= 1./24.;
//   }
//   // divided by 2 g_S^2
//   return 0.5*output/norm(FFGVertex_->norm());
				      return 1;
}

void MEPP2CharginoCharginoPowheg::constructVertex(tSubProPtr sub) {
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





