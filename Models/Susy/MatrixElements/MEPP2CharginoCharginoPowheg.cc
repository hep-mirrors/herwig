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
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
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
  MEPP2GauginoGauginoPowheg::doinit();
  // get the photon and Z ParticleData objects
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  gamma_ = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2CharginoCharginoPowheg::doinit" 
			  << Exception::abortnow;
  // do the initialisation (see Herwig::SusyBase Class)
  FFZVertex_ = hwsm->vertexFFZ();
  FFPVertex_ = hwsm->vertexFFP();
  FFGVertex_ = hwsm->vertexFFG();
  CCZVertex_ = hwsm->vertexCCZ();
  CFSVertex_ = hwsm->vertexCFSF();
  GSSVertex_ = hwsm->vertexGSFSF();
}

Selector<const ColourLines *>
MEPP2CharginoCharginoPowheg::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines c1("1 -2"), c2("1 2 -3");
  Selector<const ColourLines *> sel;
  if(abs(diag->id())==1 || abs(diag->id())==2)
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEPP2CharginoCharginoPowheg::getDiagrams() const {
  // loop over the processes we need
  tcPDPtr chi[2] = {getParticleData(1000024),getParticleData(1000037)};
  tcPDPtr chib[2];
  for(unsigned int ix=0;ix<2;++ix)
    chib[ix] = chi[ix]->CC();
  for(int i = 1; i <= maxFlavour_; ++i) {
    tcPDPtr q  = getParticleData(i);
    tcPDPtr qb = q->CC();
    tcPDPtr qL = getParticleData(1000000+i);
    tcPDPtr qR = getParticleData(2000000+i);

    for(unsigned int ix=0;ix<2;++ix){
      for(unsigned int jx=0;jx<2;++jx){
	if(process_==0 || process_== 2*ix+jx+1 ){
	  // Z-mediated s-channel
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, Z0_,
		       3, chi[ix], 3, chib[jx], -1)));
	  // photon mediated s-channel
	  add(new_ptr((Tree2toNDiagram(2), q, qb, 1, gamma_,
		       3, chi[ix], 3, chib[jx], -2)));
	  // down type
	  if(i%2==1) {
	    // ~qL mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
			 3, chi[ix], 1, chib[jx], -3)));
	    // ~qR mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
			 3, chi[ix], 1, chib[jx], -4)));
	  }
	  // up type
	  else {
	    // ~qL mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qL, qb,
			 1, chi[ix], 3, chib[jx], -3)));
	    // ~qR mediated t-channel
	    add(new_ptr((Tree2toNDiagram(3), q, qR, qb,
			 1, chi[ix], 3, chib[jx], -4)));
	  }
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
  os << FFZVertex_ << FFPVertex_ << FFGVertex_ << CCZVertex_ << CFSVertex_
     << GSSVertex_ << Z0_ << gamma_ << process_ << maxFlavour_;
}

void MEPP2CharginoCharginoPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_ >> CCZVertex_ >> CFSVertex_
     >> GSSVertex_ >> Z0_ >> gamma_ >> process_ >> maxFlavour_;
}

ClassDescription<MEPP2CharginoCharginoPowheg> 
MEPP2CharginoCharginoPowheg::initMEPP2CharginoCharginoPowheg;
// Definition of the static class description member.

void MEPP2CharginoCharginoPowheg::Init() {

  static ClassDocumentation<MEPP2CharginoCharginoPowheg> documentation
    ("MEPP2CharginoCharginoPowheg implements the ME calculation"
     " of the fermion-antifermion to chargino-chargino"
     " hard process.");

  static Switch<MEPP2CharginoCharginoPowheg,int> interfaceProcess
    ("Process",
     "Which processes to generate",
     &MEPP2CharginoCharginoPowheg::process_, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Generate all the processes"
     " (i.e. all combinations of chargino"
     "mass eigenstate pairs.)",
     0);
  static SwitchOption interfaceProcessChargino11
    (interfaceProcess,
     "chi+1chi-1",
     "Only produce chi+1, chi-1 pairs.",
     1);
  static SwitchOption interfaceProcessChargino12
    (interfaceProcess,
     "chi+1chi-2",
     "Only produce chi+1, chi-2 pairs.",
     2);
  static SwitchOption interfaceProcessChargino21
    (interfaceProcess,
     "chi+2chi-1",
     "Only produce chi+2, chi-1 pairs.",
     3);
  static SwitchOption interfaceProcessChargino22
    (interfaceProcess,
     "chi+2chi-2",
     "Only produce chi+2, chi-2 pairs.",
     4);

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
    if ( diags[i]->id() == -1)     sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 ) sel.insert(meInfo()[1], i);
    else if ( diags[i]->id() == -3 ) sel.insert(meInfo()[2], i);
    else if ( diags[i]->id() == -4 ) sel.insert(meInfo()[3], i);
  }
  return sel;
}

NLODrellYanBase::Singular MEPP2CharginoCharginoPowheg::virtualME() const {
  Singular output;
  output.eps2 = -2;
  output.eps1 = -3;
  //output.finite =-8.+sqr(Constants::pi);
  //output.finite *= loWeight();


  // average of left/right squark masses
  tcPDPtr squarkL, squarkR;
  if (abs(mePartonData()[0]->id())%2==0) {
    squarkL = getParticleData(1000000+abs(mePartonData()[0]->id())-1);
    squarkR = getParticleData(2000000+abs(mePartonData()[0]->id())-1);
  }
  else {
    squarkL = getParticleData(1000000+abs(mePartonData()[0]->id())+1);
    squarkR = getParticleData(2000000+abs(mePartonData()[0]->id())+1);
  }
  //tcPDPtr squarkL = getParticleData(1000000+mePartonData()[0]->id());
  //tcPDPtr squarkR = getParticleData(2000000+mePartonData()[0]->id());
  Energy  ms = 0.5*(squarkL->mass()+squarkR->mass());
  // boson mass
  Energy2 mz2 = sqr(Z0_->mass());
  // mandelstam variables
  Energy2 sz = sHat()-mz2;
  // I
  Complex ii(0.,1.); 
  // couplings of the vector boson
  FFVVertexPtr vertex = dynamic_ptr_cast<FFVVertexPtr>(FFPVertex_);
  vertex->setCoupling(scale(),mePartonData()[0]->CC(),
		      mePartonData()[1]->CC(),gamma_);
  double ee = vertex->electroMagneticCoupling(scale());
  double v1 = 0.5*real((vertex->left()+vertex->right())*vertex->norm())/ee;
  vertex = dynamic_ptr_cast<FFVVertexPtr>(FFZVertex_);
  vertex->setCoupling(scale(),mePartonData()[0]->CC(),
		      mePartonData()[1]->CC(),Z0_);
  ee = vertex->electroMagneticCoupling(scale());
  double v2 = 0.5*real((vertex->left()+vertex->right())*vertex->norm())/ee;
  double a2 = 0.5*real((vertex->left()-vertex->right())*vertex->norm())/ee;
  vertex = dynamic_ptr_cast<FFVVertexPtr>(CCZVertex_);
  vertex->setCoupling(scale(),mePartonData()[2],
		      mePartonData()[3],Z0_);
  ee = vertex->electroMagneticCoupling(scale());
  int v1w;
  Complex v2w;

  if(mePartonData()[2]->id() == -mePartonData()[3]->id()){
    v1w = -1.;
    v2w = -0.5*(vertex->left()+vertex->right())*vertex->norm()/ee;}
  else {
    v1w = 0.;
    v2w = 0.;
  }
  Complex a2w = 0.5*(vertex->left()-vertex->right())*vertex->norm()/ee;
  // left squark couplings
  vector<Complex> Cl(4,0.);
  FFSVertexPtr vertex2=dynamic_ptr_cast<FFSVertexPtr>(CFSVertex_);
  vertex2->setCoupling(scale(),mePartonData()[0]->CC(),mePartonData()[2],squarkL);
  ee = vertex2->electroMagneticCoupling(scale());
  Cl[0] = -0.5*vertex2->left() *vertex2->norm()/ee;
  vertex2->setCoupling(scale(),mePartonData()[1]->CC(),mePartonData()[3],squarkL->CC());
  Cl[1] =  -0.5*conj(vertex2->right()*vertex2->norm())/ee;

  tcPDPtr ccQuark;
  tcPDPtr ccSquark=getParticleData(1000000+abs(mePartonData()[0]->id()));
  if (abs(mePartonData()[0]->id())%2==0) {
    ccQuark = getParticleData(mePartonData()[0]->id()-1);
  }
  else {
    ccQuark = getParticleData(mePartonData()[0]->id()+1);
  }
  vertex2->setCoupling(scale(),ccQuark->CC(),mePartonData()[2],ccSquark);
  ee = vertex2->electroMagneticCoupling(scale());
  Cl[2] = -0.5*vertex2->left() *vertex2->norm()/ee;
  vertex2->setCoupling(scale(),ccQuark,mePartonData()[3],ccSquark->CC());
  Cl[3] =  -0.5*conj(vertex2->right()*vertex2->norm())/ee;


  if(mePartonData()[0]->id()%2!=0) {
    swap(Cl[2],Cl[0]);
    swap(Cl[3],Cl[1]);
  }

  // right squark couplings
  vector<Complex> Cr(4,0.);

  // s-channel
  vector<double> Cs(4,0.);
  Cs[0] = sqr(v1*v1w) + 2.*sHat() * real( v1*v1w*v2*v2w )/sz +
    sqr(sHat()/sz)*( sqr(v2)+sqr(a2) )*( norm(v2w) + norm(a2w) );
  Cs[1] = sqr(v1*v1w) + 2.0 * sHat() * real( v1*v1w*v2*v2w )/sz +
    sqr(sHat()/sz)*( sqr(v2)+sqr(a2) )*( norm(v2w) - norm(a2w) );
  Cs[2] = real( v1*v1w*a2*a2w ) +
    sHat()/sz*2.*a2*v2*real(a2w*conj(v2w));
  Cs[3] = 0.;
  // t-channel
  vector<Complex> Ct(4,0.);
//   Ct[0] = v1*v1w + sHat()/sz*( v2+a2 )*( v2w-a2w );
//   Ct[1] = v1*v1w + sHat()/sz*( v2+a2 )*( v2w+a2w );
//   Ct[2] = v1*v1w + sHat()/sz*( v2-a2 )*( v2w+a2w );
//   Ct[3] = v1*v1w + sHat()/sz*( v2-a2 )*( v2w-a2w );


  Ct[0] = sHat()/sz*( v2+a2 )*( v2w-a2w );
  Ct[1] = sHat()/sz*( v2+a2 )*( v2w+a2w );
  Ct[2] = sHat()/sz*( v2-a2 )*( v2w+a2w );
  Ct[3] = sHat()/sz*( v2-a2 )*( v2w-a2w );
  Ct[0] += v1*v1w;
  Ct[1] += v1*v1w;
  Ct[2] += v1*v1w;
  Ct[3] += v1*v1w;

  vector<Complex> Cv(4,0.);
  Cv[0] = Cl[2] / Cl[0];
  Cv[1] = Cl[3] / Cl[1];
  Cv[2] = Cl[0] / Cl[2];
  Cv[3] = Cl[1] / Cl[3];

  // weird rescaling factors
  for(unsigned int ix=0;ix<4;++ix) {
    Cs[ix] *= sqr(4.*Constants::pi);
    Ct[ix] *=     4.*Constants::pi ;
    Cl[ix] *= 2.0 * sqrt(Constants::pi);
    Cr[ix] *= 2.0 * sqrt(Constants::pi);
  }
  // finite piece
  output.finite = finiteVirtual(ms,mz2,Cl,Cr,Cs,Ct,Cv);


  return output;
}

double MEPP2CharginoCharginoPowheg::
qqbarME(vector<SpinorWaveFunction>    & sp ,
	vector<SpinorBarWaveFunction> & sbar ,
	vector<SpinorWaveFunction>    & spout ,
	vector<SpinorBarWaveFunction> & sbarout ,
	//ScalarWaveFunction & interqL,ScalarWaveFunction & interqR,
	bool first) const {
  // scale for the process
  const Energy2 q2(scale());
  tcPDPtr squark[2];
  if (abs(mePartonData()[0]->id())%2==0) {
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())-1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())-1);
  }
  else {
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())+1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())+1);
  }
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
  vector<double> me(4, 0.);
  double me2(0.);
  // storage of the individual diagrams
  vector<Complex> diag(4, 0.);
  ProductionMatrixElement pme(PDT::Spin1Half, PDT::Spin1Half, 
			      PDT::Spin1Half, PDT::Spin1Half);
  unsigned int propopt;
  if(status()==RealQG || status()==RealQBarG){
    propopt = 7;
  }
  else propopt = 3;
  // loop over the helicities and calculate the matrix elements
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      VectorWaveFunction interZ = FFZVertex_->
	evaluate(q2, 1, Z0_, sp[if1],sbar[if2]);
      VectorWaveFunction interP = FFPVertex_->
	evaluate(q2, 1, gamma_, sp[if1],sbar[if2]);
      for(unsigned int of1 = 0; of1 < 2; ++of1) {
	for(unsigned int of2 = 0; of2 < 2; ++of2) {
	  // s-channel
  	  diag[0] = CCZVertex_->evaluate(q2, spout[of1],  sbarout[of2], interZ);
 	  if(spout[of1].particle()->id()==-sbarout[of2].particle()->id())
 	    diag[1] = CCZVertex_->evaluate(q2, spout[of1],  sbarout[of2], interP);
	  // t-channel squark exchanges	  
	  for(unsigned int iq=0;iq<2;++iq) {
	    if(abs(mePartonData()[0]->id())%2==0) {
	      ScalarWaveFunction intersq = CFSVertex_->
		evaluate(q2, propopt, squark[iq], sp[if1], spoutconj[of1]);
	      diag[iq+2] = 
		-CFSVertex_->evaluate(q2, sbaroutconj[of2], sbar[if2], intersq);
	    }
	    else {
	      ScalarWaveFunction intersq = CFSVertex_->
		evaluate(q2, propopt, squark[iq], sp[if1], sbarout[of2]);
	      diag[iq+2] = 
		CFSVertex_->evaluate(q2, spout[of1], sbar[if2], intersq);
	    }
	  }
	  // individual diagrams
	  for(unsigned int id=0;id<4;++id){
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
    DVector save(4);
    for(DVector::size_type ix = 0; ix < 4; ++ix)
      save[ix] = me[ix]/12.;
    meInfo(save);
    me_.reset(pme);
  }
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
  // outgoing chargino wavefunctions
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
  vector<SpinorWaveFunction> sp(2);
  vector<SpinorBarWaveFunction> sbar(2);
  vector<VectorWaveFunction> gluon(2);
  // squarks for the t-channel
  tcPDPtr squark[2];
  if (abs(mePartonData()[0]->id())%2==0) {
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())-1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())-1);
  }
  else {
    squark[0] = getParticleData(1000000+abs(mePartonData()[0]->id())+1);
    squark[1] = getParticleData(2000000+abs(mePartonData()[0]->id())+1);
  }
  // wavefunctions for the q qbar -> chi chi g process
  if(particles[0]->id()==-particles[1]->id()) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[1],particles[1],  i,incoming);
      gluon[i]= VectorWaveFunction   (momenta[4],particles[4],2*i,outgoing);
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
  // g q
  else if(particles[0]->id()>0 &&
	  particles[1]->id()==ParticleID::g) {
    for( unsigned int i = 0; i < 2; ++i ) {
      sp[i]   = SpinorWaveFunction   (momenta[0],particles[0],  i,incoming);
      sbar[i] = SpinorBarWaveFunction(momenta[4],particles[4],  i,outgoing);
      gluon[i]= VectorWaveFunction   (momenta[1],particles[1],2*i,incoming);
    }
  }
  else {
    for(unsigned int ix=0;ix<particles.size();++ix) {
      cerr << particles[ix]->PDGName() << " " << momenta[ix]/GeV << "\n";
    }
    assert(false);
  }
  // wavefunctions for the outgoing charginos
  vector<SpinorWaveFunction> spout(2),sbaroutconj(2);
  vector<SpinorBarWaveFunction> sbarout(2),spoutconj(2);
  for( unsigned int i = 0; i < 2; ++i ) {
    spout[i]       =    SpinorWaveFunction(momenta[2], particles[2], i,
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
  vector<Complex> diag(10,0.);
  Energy2 q2 = scale();
  unsigned int propopt;
  if(status()==RealQG || status()==RealQBarG){
    propopt = 7;
  }
  else propopt = 3;
  for(unsigned int ihel1=0;ihel1<2;++ihel1) {
    for(unsigned int ihel2=0;ihel2<2;++ihel2) {
      for(unsigned int ohel1=0;ohel1<2;++ohel1) {
	SpinorWaveFunction    inters = FFGVertex_->evaluate(q2, 5, sp[ihel1].particle()->CC(),
							    sp[ihel1], gluon[ohel1]);
	SpinorBarWaveFunction interb = FFGVertex_->evaluate(q2,5,sbar[ihel2].particle()->CC(),
							    sbar[ihel2],gluon[ohel1]);
	VectorWaveFunction interV[4] = 
	  {FFZVertex_->evaluate(q2, 1, Z0_, inters,sbar[ihel2]),
	   FFZVertex_->evaluate(q2, 1, Z0_, sp[ihel1],interb),
	   FFPVertex_->evaluate(q2, 1, gamma_, inters,sbar[ihel2]),
	   FFPVertex_->evaluate(q2, 1, gamma_, sp[ihel1],interb)};
	for(unsigned int ohel2=0;ohel2<2;++ohel2) {
	  for(unsigned int ohel3=0;ohel3<2;++ohel3) {
	    // s-channel diagrams
	    // first Z diagram (emission from quark)
	    diag[0] = CCZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[0]);
	    // second Z diagram (emission from anti-quark)
	    diag[1] = CCZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[1]);
	    //if(mePartonData[2]()->id() == -mePartonData[3]()->id()){
	    if(spout[ohel1].particle()->id()==-sbarout[ohel2].particle()->id()){
	      // first photon diagram (emission from quark)
	      diag[2] = CCZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[2]);
	      // second photon diagram (emission from anti-quark)
	      diag[3] = CCZVertex_->evaluate(q2, spout[ohel2],  sbarout[ohel3], interV[3]);
	    }
	    // t-channel squark exchanges
	    ScalarWaveFunction intersq,intersq2;	  
	    for(unsigned int iq=0;iq<2;++iq) {
	      if(   (2.*squark[iq]->mass())/(momenta[2].mass()+momenta[3].mass()) > 20.){
		diag[3*iq+4]=0.;
		diag[3*iq+5]=0.;
		diag[3*iq+6]=0.;
	      }
	      else{
		// u-type
		if(abs(mePartonData()[0]->id())%2==0) {
		  // emission from quark
		  intersq = CFSVertex_->
		    evaluate(q2, propopt, squark[iq], inters, spoutconj[ohel2]);
		  diag[3*iq+4] = 
		    -CFSVertex_->evaluate(q2, sbaroutconj[ohel3], sbar[ihel2], intersq);
		  // emission antiquark
		  intersq = CFSVertex_->
		    evaluate(q2, propopt, squark[iq], sp[ihel1], spoutconj[ohel2]);
		  diag[3*iq+5] = 
		    -CFSVertex_->evaluate(q2, sbaroutconj[ohel3], interb, intersq);
		  // emission from intermediate
		  intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq],gluon[ohel1],intersq);
		  diag[3*iq+6] = 
		    -CFSVertex_->evaluate(q2, sbaroutconj[ohel3], sbar[ihel2], intersq2);
		}
		// down type
		else {
		  // emission quark
		  intersq = CFSVertex_->
		    evaluate(q2, propopt, squark[iq], inters, sbarout[ohel3]);
		  diag[3*iq+4] = 
		    CFSVertex_->evaluate(q2, spout[ohel2], sbar[ihel2], intersq);
		  // emission antiquark
		  intersq = CFSVertex_->
		    evaluate(q2, propopt, squark[iq], sp[ihel1], sbarout[ohel3]);
		  diag[3*iq+5] = 
		    CFSVertex_->evaluate(q2, spout[ohel2], interb, intersq);
		  // emission from intermediate
		  intersq2 = GSSVertex_->evaluate(q2,propopt,squark[iq],gluon[ohel1],intersq);
		  diag[3*iq+6] = 
		    CFSVertex_->evaluate(q2, spout[ohel2], sbar[ihel2], intersq2);
		}
	      }
	    }
	    // add them up
	    Complex total = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    output += norm(total);
	  }
	}
      }
    }
  }


  // strong coupling
  double gs2 = norm(FFGVertex_->norm());
  
  
  double totcount = 0.;
//   Energy sqmom,sqmass;
//   bool first = false;
//   bool second = false;


  // subtract the on-shell squark decay if neccessary
  if(abs(particles[4]->id())<=6) {
    Energy roots = (momenta[0]+momenta[1]).m();
    // off-shell masses of the squarks
    Energy2 msq1 = (momenta[2]+momenta[4]).m2();
    Energy2 msq2 = (momenta[3]+momenta[4]).m2();
    Lorentz5Momentum pgluon = 
      particles[0]->id() == ParticleID::g ? momenta[0] : momenta[1];
    tcPDPtr quark = particles[0]->id() == ParticleID::g ?
      particles[1] : particles[0];
    FFSVertexPtr vertex = dynamic_ptr_cast<FFSVertexPtr>(CFSVertex_);
    //FFSVertexPtr vertexd = dynamic_ptr_cast<FFSVertexPtr>(CFSVertex_);
    double Cf = 4./3.;
    double Nc = 3.;
    Energy2 sh = (momenta[0]+momenta[1]).m2();
    for(unsigned int ix=0;ix<2;++ix) {
      Energy2 sqwidth2 = sqr(squark[ix]->width());

      // is on-shell mass allowed for first chargino and squark
      if(  (roots >= squark[ix]->mass() + momenta[2].mass() )
	   && (squark[ix]->mass() > momenta[3].mass())
	   && ((quark->iCharge() > 0. && particles[2]->iCharge() > 0.) || 
	       (quark->iCharge() < 0. && particles[2]->iCharge() < 0.))){

// 	first = true;

// 	sqmom = (momenta[3]+momenta[4]).m();
// 	sqmass = squark[ix]->mass();


	// Create a counterterm for the squark decay pole.
	Energy2 mcharg2 = sqr(momenta[2].mass()); 
	Energy2 t3 = (pgluon-momenta[3]-momenta[4]).m2()-msq2;
	Energy2 u4 = (pgluon-momenta[2]).m2()-mcharg2;
	vertex->setCoupling(q2,quark,particles[2],squark[ix]);
 	double a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	double sqprod2 = (1./8.) * gs2 * Cf * Nc * a2 *
	  (-u4/sh - (2*(msq2-mcharg2)*u4)/sh/t3 * (1+mcharg2/u4+msq2/t3));
	//vertex->setCoupling(q2,quark,particles[3],squark[ix]);
	vertex->setCoupling(q2,particles[4],particles[3],squark[ix]);

// 	if ((quark->iCharge() > 0. && particles[2]->iCharge() < 0.) || 
// 	    (quark->iCharge() < 0. && particles[2]->iCharge() > 0.)){
// 	  cout << "1st CT:" << endl;
// 	  cout << "Incoming partons: " << particles[0]->id() << " " << particles[1]->id() << endl;
// 	  cout << "Outgoing Charginos: " << particles[2]->id() << "\t" << particles[3]->id() << endl;
// 	  cout << "Born Chargino: " << particles[2]->id() << " " << "Decay Chargino " << particles[3]->id() << "\n" << endl;}

 	a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	Energy2 sqdecay2 = 4. * a2 * (msq2-sqr(particles[3]->mass()));
	Energy4 denom = sqr(msq2-sqr(squark[ix]->mass())) + 
	  sqr(squark[ix]->mass())*sqwidth2;
	double sqcounter = sqprod2 * sqdecay2 * UnitRemoval::E2 / denom;

	if( ((2.*squark[ix]->mass())/(momenta[2].mass()+momenta[3].mass()) > 100.)
	    && (squark[ix]->mass()>1.e4*GeV))
	  sqcounter = 0;

	totcount += sqcounter;
      }

      // is on-shell mass allowed for second chargino and squark
      if((roots >= squark[ix]->mass() + momenta[3].mass() )
	 && (squark[ix]->mass() > momenta[2].mass())
	 && ((quark->iCharge() > 0. && particles[3]->iCharge() > 0.) || 
	     (quark->iCharge() < 0. && particles[3]->iCharge() < 0.))){
	
	
// 	second = true;

// 	sqmom = (momenta[2]+momenta[4]).m();
// 	sqmass = squark[ix]->mass();


	Energy2 mcharg2 = sqr(momenta[3].mass());
	Energy2 t3 = (pgluon-momenta[2]-momenta[4]).m2()-msq1;
	Energy2 u4 = (pgluon-momenta[3]).m2()-mcharg2;
	vertex->setCoupling(q2,quark,particles[3],squark[ix]);
 	double a2 = norm(vertex->norm())*
 	  (norm(vertex->left())+norm(vertex->right()));
	double sqprod2 = (1./8.) * gs2 * Cf * Nc * a2 *
	  (-u4/sh - (2*(msq1-mcharg2)*u4)/sh/t3 * (1+mcharg2/u4+msq1/t3));
	//vertex->setCoupling(q2,quark,particles[2],squark[ix]);
	vertex->setCoupling(q2,particles[4],particles[2],squark[ix]);

// 	if ((quark->iCharge() > 0. && particles[3]->iCharge() < 0.) || 
// 	    (quark->iCharge() < 0. && particles[3]->iCharge() > 0.)){
// 	cout << "2nd CT:" << endl;
// 	cout << "Incoming partons: " << particles[0]->id() << " " << particles[1]->id() << endl;
// 	cout << "Outgoing Charginos: " << particles[2]->id() << "\t" << particles[3]->id() << endl;
// 	cout << "Born Chargino: " << particles[3]->id() << " " << "Decay Chargino: " << particles[2]->id() << "\n" << endl;}

	a2 = norm(vertex->norm())*
	  (norm(vertex->left())+norm(vertex->right()));
	Energy2 sqdecay2 = 4. * a2 * (msq1-sqr(particles[2]->mass()));
	Energy4 denom = sqr(msq1-sqr(squark[ix]->mass())) + 
	  sqr(squark[ix]->mass())*sqwidth2;
	double sqcounter = sqprod2 * sqdecay2 * UnitRemoval::E2 / denom;

	if( ((2.*squark[ix]->mass())/(momenta[2].mass()+momenta[3].mass()) > 100.)
	    && (squark[ix]->mass()>1.e4*GeV))
	  sqcounter = 0;

	totcount += sqcounter;
      }
    }
  }



//   if(abs(1.-abs(sqmom/sqmass)) < 1.e-4){
//     cout << "Counterterms on:\t" << first << "\t" << second << endl;
//     cout << particles[0]->id() << "\t" << particles[1]->id() << endl;
//     cout << particles[2]->id() << "\t" << particles[3]->id() << "\t" << particles[4]->id() << endl;
//     cout << sqmom/GeV << "\t" << sqmass/GeV << "\t" << "\t ratio: " << output << "\t" << totcount << "\t" << output/totcount << "\n" << endl;
//  }
  

  output -= totcount;


  // colour and spin factors
  if(particles[0]->id()==-particles[1]->id()) {
    output *= 1./9.;
  }
  else  {
    output *= 1./24.;
  }
  // divided by 2 g_S^2
  return 0.5*output/norm(FFGVertex_->norm());
}

void MEPP2CharginoCharginoPowheg::constructVertex(tSubProPtr ) {
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





