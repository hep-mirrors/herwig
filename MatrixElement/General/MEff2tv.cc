// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEff2tv class.
//

#include "MEff2tv.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MEff2tv::clone() const {
  return new_ptr(*this);
}

IBPtr MEff2tv::fullclone() const {
  return new_ptr(*this);
}

void MEff2tv::persistentOutput(PersistentOStream & os) const {
  os << fermion_ << vector_ << fourPoint_;
}

void MEff2tv::persistentInput(PersistentIStream & is, int) {
  is >> fermion_ >> vector_ >> fourPoint_;
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin2    , PDT::Spin1);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEff2tv,GeneralHardME>
describeHerwigMEff2tv("Herwig::MEff2tv", "Herwig.so");

void MEff2tv::Init() {

  static ClassDocumentation<MEff2tv> documentation
    ("The MEff2tv class implements the general matrix element for "
     "fermion-antifermion -> tensor vector");

}

double MEff2tv::me2() const {
  // first setup  wavefunctions for external particles
  SpinorVector sp(2);
  SpinorBarVector sbar(2);
  VBVector vec(3);
  TBVector ten(5);
  bool tMass = meMomenta()[2].mass()!=ZERO;
  bool vMass = meMomenta()[3].mass()!=ZERO;
  for( unsigned int i = 0; i < 5; ++i ) {
    if(i<2) {
      sp[i]   = SpinorWaveFunction   (rescaledMomenta()[0], 
				      mePartonData()[0], i, incoming);
      sbar[i] = SpinorBarWaveFunction(rescaledMomenta()[1], 
				      mePartonData()[1], i, incoming);
    }
    if( tMass || i==0 || i==4) {
      ten[i] = TensorWaveFunction(rescaledMomenta()[2], mePartonData()[2],i , 
				  outgoing);
    }
    if(i<3 && (i!=1||vMass) ) {
      vec[i] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3],i , 
 				  outgoing);
    }
  }
  // calculate the ME
  double full_me(0.);
  ffb2tvHeME(sp, sbar, ten, vec, full_me,true);  
  // debugging tests if needed
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  // return the answer
  return full_me;
}

void MEff2tv::doinit() {
  GeneralHardME::doinit();
  fermion_   .resize(numberOfDiags());
  vector_    .resize(numberOfDiags());
  fourPoint_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1Half, PDT::Spin1Half, 
			   PDT::Spin2    , PDT::Spin1);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1Half)
	throw InitException() << "MEff2tv:doinit() - Cannot find correct "
			      << "t-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      if( current.ordered.second ) 
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.first), 
		    dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.second));
      else
	fermion_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractFFTVertexPtr>(current.vertices.second), 
		    dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first));
    }
    else if(current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1)
	throw InitException() << "MEff2tv:doinit() - Cannot find correct "
			      << "s-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      vector_[i] = 
	make_pair(dynamic_ptr_cast<AbstractFFVVertexPtr>(current.vertices.first), 
		  dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.second));
    }
    else if(current.channelType == HPDiagram::fourPoint) {
      fourPoint_[i] = 
	dynamic_ptr_cast<AbstractFFVTVertexPtr>(current.vertices.first);
    }
  }
}

ProductionMatrixElement MEff2tv::
ffb2tvHeME(SpinorVector & sp, SpinorBarVector & sb,
	   TBVector & ten, VBVector & vec, 
	   double & me2,bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool tMass = meMomenta()[2].mass() != ZERO;
  bool vMass = meMomenta()[3].mass() != ZERO;
  // flow over the helicities and diagrams
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int it=0; it<5; ++it) {
	if( (it>0&&it<4) && !tMass ) continue;
	for(unsigned int iv=0; iv<3;++iv) {
	  if(iv==1&&!vMass) continue;
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr internal(current.intermediate);	
	    if(current.channelType == HPDiagram::tChannel) {
	      if(current.ordered.second) {
		if(internal->CC()) internal = internal->CC();
		SpinorBarWaveFunction interFB = fermion_[ix].second->
		  evaluate(q2,5,internal,sb[if2],vec[iv]);
		diag = fermion_[ix].first->
		  evaluate(q2,sp[if1],interFB,ten[it]);
	      }
	      else {
		SpinorWaveFunction interF = fermion_[ix].second->
		  evaluate(q2,5,internal,sp[if1],vec[iv]);
		diag = fermion_[ix].first->
		  evaluate(q2,interF,sb[if2],ten[it]);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      VectorWaveFunction interV = vector_[ix].first->
		evaluate(q2, 1, internal, sp[if1], sb[if2],vec[iv].mass());
	      diag = vector_[ix].second->evaluate(q2, interV, vec[iv],ten[it],
						  vec[iv].mass());
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      diag = fourPoint_[ix]->
		evaluate(q2,sp[if1],sb[if2],vec[iv],ten[it]);
	    }
	    // diagram
	    me[ix] += norm(diag);
	    diagramME()[ix](if1,if2,it,iv) = diag;
	    // contributions to the different colour flows
	    for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](if1,if2,it,iv) = flows[iy];
	  // contribution to the squared matrix element
	  for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
	    for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
	      me2 += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
	  // contribution to the colour flow
	  for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) {
	    flow[ii] += getColourFactors()[ii][ii]*norm(flows[ii]);
	  }
	}
      }
    }
  }
  // if not computing the cross section return the selected colour flow
  if(!first) return flowME()[colourFlow()];
  me2 = selectColourFlow(flow,me,me2);
  return flowME()[colourFlow()];
}

void MEff2tv::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id4 = mePartonData()[3]->id();
  if(id1==-id2&&id1<=5&&id4==ParticleID::g) {
    unsigned int iloc(0);
    for(;iloc<vector_.size();++iloc)
      if(vector_[iloc].first) break;
    double gs = abs(vector_[iloc].first->norm());
    InvEnergy kappa = abs(vector_[iloc].second->norm())*UnitRemoval::InvE;
    Energy2 mg2 = sqr(meMomenta()[2].mass());
    double anal = sqr(gs)*sqr(kappa)/36.*(4.*uHat()*tHat()+sHat()*mg2)*
      (sqr(tHat()-mg2)+sqr(uHat()-mg2))/sHat()/tHat()/uHat();
    double diff = abs((anal - me2)/(anal+me2));
    if( diff > 1e-4 ) {
      generator()->log() 
	<< mePartonData()[0]->PDGName() << "," 	
	<< mePartonData()[1]->PDGName() << "->"
	<< mePartonData()[2]->PDGName() << ","
	<< mePartonData()[3]->PDGName() << "   difference: " 
	<< setprecision(10) << diff << "  ratio: " << anal/me2 
	<< '\n';
    }
  }
  else if(id1==-id2&&id1==ParticleID::eminus&&id4==ParticleID::gamma) {
    unsigned int iloc(0);
    for(;iloc<vector_.size();++iloc)
      if(vector_[iloc].first) break;
    double gs = abs(vector_[iloc].first->norm());
    InvEnergy kappa = abs(vector_[iloc].second->norm())*UnitRemoval::InvE;
    Energy2 mg2 = sqr(meMomenta()[2].mass());
    double anal = sqr(gs)*sqr(kappa)/16./tHat()/uHat()/sHat()*
      (4.*uHat()*tHat()+mg2*sHat())*(sqr(uHat()-mg2)+sqr(tHat()-mg2));
    double diff = abs((anal - me2)/(anal+me2));
    if( diff > 1e-4 ) {
      generator()->log() 
	<< mePartonData()[0]->PDGName() << "," 	
	<< mePartonData()[1]->PDGName() << "->"
	<< mePartonData()[2]->PDGName() << ","
	<< mePartonData()[3]->PDGName() << "   difference: " 
	<< setprecision(10) << diff << "  ratio: " << anal/me2 
	<< '\n';
    }    
  }
  else if(id1==-id2&&id1==ParticleID::eminus&&id4==ParticleID::Z0) {
    unsigned int iloc(0);
    for(;iloc<vector_.size();++iloc)
      if(vector_[iloc].first) break;
    double gs = abs(vector_[iloc].first->norm());
    InvEnergy kappa = abs(vector_[iloc].second->norm())*UnitRemoval::InvE;
    Energy2 mg2 = sqr(meMomenta()[2].mass());
    Energy2 mz2 = sqr(meMomenta()[3].mass());
    double sw2 = SM().sin2ThetaW();
    double anal = sqr(gs)*sqr(kappa)/48./16./sw2/(1.-sw2)*
      2.*(1.-4.*sw2+8.*sqr(sw2))/sqr(tHat())/sqr(uHat())/sqr(sHat()-mz2)*
      (8.*pow<3,1>(mz2)*uHat()*tHat()*(3.*mg2*(mg2-uHat()-tHat())+4.*uHat()*tHat())
       +2.*sqr(mz2)*uHat()*tHat()*(27.*pow<3,1>(mg2)-42.*sqr(mg2)*(uHat()+tHat())
				   +15.*mg2*(sqr(uHat())+sqr(tHat()))
				   +80.*mg2*uHat()*tHat()
				   -28.*(sqr(uHat())*tHat()+uHat()*sqr(tHat())))
       +mz2*(3.*pow<4,1>(mg2)*(-sqr(uHat())-sqr(tHat())+12.*uHat()*tHat())
	     +6.*pow<3,1>(mg2)*(pow<3,1>(uHat())
				-12.*(sqr(uHat())*tHat()+uHat()*sqr(tHat()))
				+pow<3,1>(tHat()))
	     +3.*sqr(mg2)*(-pow<4,1>(uHat())+14.*pow<3,1>(uHat())*tHat()
			   +62.*sqr(uHat()*tHat())+14.*uHat()*pow<3,1>(tHat())
			   -pow<4,1>(tHat()))
	     +6.*mg2*(-pow<4,1>(uHat())*tHat()
		      -23.*(pow<3,1>(uHat())*sqr(tHat())+pow<3,1>(tHat())*sqr(uHat()))
		      -uHat()*pow<4,1>(tHat()))
	     +36.*(pow<4,1>(uHat())*sqr(tHat())+pow<4,1>(tHat())*sqr(uHat()))
	     +52.*pow<3,1>(tHat()*uHat()))
       +3.*tHat()*uHat()*(-mg2+uHat()+tHat())*
       (-sqr(mg2)+mg2*(uHat()+tHat())-4.*uHat()*tHat())*
       (2.*sqr(mg2)-2.*mg2*(uHat()+tHat())+sqr(uHat())+sqr(tHat())));
    double diff = abs((anal - me2)/(anal+me2));
    if( diff > 1e-4 ) {
      generator()->log() 
	<< mePartonData()[0]->PDGName() << "," 	
	<< mePartonData()[1]->PDGName() << "->"
	<< mePartonData()[2]->PDGName() << ","
	<< mePartonData()[3]->PDGName() << "   difference: " 
	<< setprecision(10) << diff 
	<< " anal : " << anal
	<< " code : " << me2
	<< " ratio: " << anal/me2 
	<< '\n';
    }    
  }
}

void MEff2tv::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  vector<SpinorWaveFunction> sp;
  SpinorWaveFunction(sp, ext[0], incoming, false);
  vector<SpinorBarWaveFunction> sbar;
  SpinorBarWaveFunction(sbar, ext[1], incoming, false);
  vector<VectorWaveFunction> v2;
  vector<TensorWaveFunction> t1;
  bool mc  = !(ext[2]->momentum().mass() > ZERO);
  bool md  = !(ext[3]->data()    .mass() > ZERO);
  TensorWaveFunction(t1, ext[2], outgoing, true, mc);
  VectorWaveFunction(v2, ext[3], outgoing, true, md);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  SpinorWaveFunction spr   (rescaledMomenta()[0],
			    ext[0]->dataPtr(), incoming);
  SpinorBarWaveFunction sbr(rescaledMomenta()[1],
			    ext[1]->dataPtr(), incoming);
  TensorWaveFunction tr1   (rescaledMomenta()[2],
			    ext[2]->dataPtr(), outgoing);
  VectorWaveFunction vr2   (rescaledMomenta()[3],
			    ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    spr.reset(ihel);
    sp[ihel] = spr;
    sbr.reset(ihel);
    sbar[ihel] = sbr;
    tr1.reset(4*ihel);
    t1[4*ihel] = tr1;
    vr2.reset(2*ihel);
    v2[2*ihel] = vr2;
  }
  if( !mc ) {
    for(unsigned int ihel=1;ihel<4;++ihel) {
      tr1.reset(ihel);
      t1[ihel] = tr1;
    }
  }
  if( !md ) {
    vr2.reset(1);
    v2[1] = vr2;
  }
  double dummy(0.);
  ProductionMatrixElement pme = ffb2tvHeME(sp, sbar, t1, v2,dummy,false);
  createVertex(pme,ext);
}
