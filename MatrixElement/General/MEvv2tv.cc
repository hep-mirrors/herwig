// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEvv2tv class.
//

#include "MEvv2tv.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

IBPtr MEvv2tv::clone() const {
  return new_ptr(*this);
}

IBPtr MEvv2tv::fullclone() const {
  return new_ptr(*this);
}

void MEvv2tv::persistentOutput(PersistentOStream & os) const {
  os << vector_ << fourPoint_;
}

void MEvv2tv::persistentInput(PersistentIStream & is, int) {
  is >> vector_ >> fourPoint_;
  initializeMatrixElements(PDT::Spin1, PDT::Spin1, 
			   PDT::Spin2, PDT::Spin1);
}

ClassDescription<MEvv2tv> MEvv2tv::initMEvv2tv;
// Definition of the static class description member.

void MEvv2tv::Init() {

  static ClassDocumentation<MEvv2tv> documentation
    ("The MEvv2tv class implements the general matrix element "
     "for vector vector -> tensor vector.");

}

double MEvv2tv::me2() const {
  // check tensor and outgoing vector mass
  bool tMass = meMomenta()[2].mass()!=ZERO;
  bool vMass = meMomenta()[3].mass()!=ZERO;
  // first setup  wavefunctions for external particles
  VBVector vec1(3),vec2(3),vec3(3);
  TBVector ten(5);
  for( unsigned int i = 0; i < 5; ++i ) {
    if(i!=1 &&i<3) {
      vec1[i] = VectorWaveFunction(rescaledMomenta()[0], mePartonData()[0],i , 
				   incoming);
      vec2[i] = VectorWaveFunction(rescaledMomenta()[1], mePartonData()[1],i , 
				   incoming);
    }
    if( tMass || i==0 || i==4) {
      ten[i] = TensorWaveFunction(rescaledMomenta()[2], mePartonData()[2],i , 
				  outgoing);
    }
    if(i<3 && (i!=1||vMass) ) {
      vec3[i] = VectorWaveFunction(rescaledMomenta()[3], mePartonData()[3],i , 
				  outgoing);
    }
  }
  // calculate the ME
  double full_me(0.);
  vv2tvHeME(vec1,vec2,ten,vec3,full_me,true);
  // debugging tests if needed
#ifndef NDEBUG
  if( debugME() ) debug(full_me);
#endif
  // return the answer
  return full_me;
}

void MEvv2tv::doinit() {
  GeneralHardME::doinit();
  vector_    .resize(numberOfDiags());
  fourPoint_ .resize(numberOfDiags());
  initializeMatrixElements(PDT::Spin1, PDT::Spin1, 
			   PDT::Spin2, PDT::Spin1);
  for(HPCount i = 0; i < numberOfDiags(); ++i) {
    const HPDiagram & current = getProcessInfo()[i];
    if(current.channelType == HPDiagram::tChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1)
	throw InitException() << "MEvv2tv:doinit() - Cannot find correct "
			      << "t-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
      if( current.ordered.second ) 
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractVVVVertexPtr>(current.vertices.second), 
		    dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.first ));
      else
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractVVVVertexPtr>(current.vertices.first ), 
		    dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.second));
    }
    else if (current.channelType == HPDiagram::sChannel) {
      if(current.intermediate->iSpin() != PDT::Spin1)
	throw InitException() << "MEvv2tv:doinit() - Cannot find correct "
			      << "s-channel from diagram. Vertex not cast! "
			      << Exception::runerror;
	vector_[i] = 
	  make_pair(dynamic_ptr_cast<AbstractVVVVertexPtr>(current.vertices.first ), 
		    dynamic_ptr_cast<AbstractVVTVertexPtr>(current.vertices.second));
      
    }
    else if(current.channelType == HPDiagram::fourPoint) {
      fourPoint_[i] = 
	dynamic_ptr_cast<AbstractVVVTVertexPtr>(current.vertices.first);
    }
  }
}

ProductionMatrixElement MEvv2tv::vv2tvHeME(const VBVector & vec1,
					   const VBVector & vec2,
					   const TBVector & ten,
					   const VBVector & vec3,
					   double & full_me, bool first) const {
  // scale
  const Energy2 q2(scale());
  // weights for the selection of the diagram
  vector<double> me(numberOfDiags(), 0.);
  // weights for the selection of the colour flow
  vector<double> flow(numberOfFlows(),0.);
  bool tMass = meMomenta()[2].mass() != ZERO;
  bool vMass = meMomenta()[3].mass() != ZERO;
  // flow over the helicities and diagrams
  for(unsigned int iv1 = 0; iv1 < 3; iv1+=2) {
    for(unsigned int iv2 = 0; iv2 < 3; iv2+=2) {
      for(unsigned int it=0; it<5; ++it) {
	if( (it>0&&it<4) && !tMass ) continue;
	for(unsigned int iv3=0; iv3<3;++iv3) {
	  if(iv3==1&&!vMass) continue;
	  vector<Complex> flows(numberOfFlows(),0.);
	  for(HPCount ix = 0; ix < numberOfDiags(); ++ix) {
	    Complex diag(0.);
	    const HPDiagram & current = getProcessInfo()[ix];
	    tcPDPtr internal(current.intermediate);
	    if(current.channelType == HPDiagram::tChannel) {
	      if(current.ordered.second) {
		VectorWaveFunction interV =  vector_[ix].first->
		  evaluate(q2, 1, internal, vec2[iv2], vec3[iv3]);
		diag = vector_[ix].second->
		  evaluate(q2,interV,vec1[iv1],ten[it]);
	      }
	      else {
		VectorWaveFunction interV =  vector_[ix].first->
		  evaluate(q2, 1, internal, vec3[iv3], vec1[iv1]);
		diag = vector_[ix].second->
		  evaluate(q2,interV,vec2[iv2],ten[it]);
	      }
	    }
	    else if(current.channelType == HPDiagram::sChannel) {
	      VectorWaveFunction interV = vector_[ix].first->
		evaluate(q2, 1, internal, vec1[iv1], vec2[iv2]);
	      diag = vector_[ix].second->evaluate(q2, interV, vec3[iv3],ten[it]);
	    }
	    else if(current.channelType == HPDiagram::fourPoint) {
	      diag = fourPoint_[ix]->
		evaluate(q2,vec1[iv1],vec2[iv2],vec3[iv3],ten[it]);
	    }
	    // diagram
	    me[ix] += norm(diag);
	    diagramME()[ix](iv1,iv2,it,iv3) = diag;
	    // contributions to the different colour flows
	    for(unsigned int iy = 0; iy < current.colourFlow.size(); ++iy) {
	      assert(current.colourFlow[iy].first<flows.size());
	      flows[current.colourFlow[iy].first] += 
		current.colourFlow[iy].second * diag;
	    }
	  }
	  // MEs for the different colour flows
	  for(unsigned int iy = 0; iy < numberOfFlows(); ++iy) 
	    flowME()[iy](iv1,iv2,it,iv3) = flows[iy];
	  // contribution to the squared matrix element
	  for(unsigned int ii = 0; ii < numberOfFlows(); ++ii) 
	    for(unsigned int ij = 0; ij < numberOfFlows(); ++ij)
	      full_me += getColourFactors()[ii][ij]*(flows[ii]*conj(flows[ij])).real();
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
  full_me = selectColourFlow(flow,me,full_me);
  return flowME()[colourFlow()];
}

void MEvv2tv::debug(double me2) const {
  if( !generator()->logfile().is_open() ) return;
  long id1 = mePartonData()[0]->id();
  long id2 = mePartonData()[1]->id();
  long id4 = mePartonData()[3]->id();
  if(id1 != ParticleID::g || id2 != ParticleID::g || 
     id4 != ParticleID::g ) return;
  unsigned int iloc(0);
  for(;iloc<vector_.size();++iloc)
    if(vector_[iloc].first) break;
  double gs = abs(vector_[iloc].second->norm());
  InvEnergy kappa = abs(vector_[iloc].first->norm())*UnitRemoval::InvE;
  Energy2 mg2 = sqr(meMomenta()[2].mass());
  double anal2 = 3./32.*sqr(gs)*sqr(kappa)/sHat()/tHat()/uHat()
    *(pow<4,1>(sHat()-mg2)+pow<4,1>(tHat()-mg2)+pow<4,1>(uHat()-mg2));
  double diff = abs((anal2 - me2)/(anal2+me2));
  if( diff > 1e-4 ) {
    generator()->log() 
      << mePartonData()[0]->PDGName() << "," 	
      << mePartonData()[1]->PDGName() << "->"
      << mePartonData()[2]->PDGName() << ","
      << mePartonData()[3]->PDGName() << "   difference: " 
      << setprecision(10) << diff << "  ratio: " << anal2/me2 
      << '\n';
  }
}

void MEvv2tv::constructVertex(tSubProPtr sub) {
  ParticleVector ext = hardParticles(sub);
  // set wave functions with real momenta
  VBVector v1, v2;
  VectorWaveFunction(v1, ext[0], incoming, false, true);
  VectorWaveFunction(v2, ext[1], incoming, false, true);
  //function to calculate me2 expects massless incoming vectors
  // and this constructor sets the '1' polarisation at element [2] 
  //in the vector
  bool mc  = !(ext[2]->momentum().mass() > ZERO);
  bool md  = !(ext[3]->data()    .mass() > ZERO);
  VBVector v4;
  vector<TensorWaveFunction> t3;
  TensorWaveFunction(t3, ext[2], outgoing, true, mc);
  VectorWaveFunction(v4, ext[3], outgoing, true, md);
  // Need to use rescale momenta to calculate matrix element
  setRescaledMomenta(ext);
  // wave functions with rescaled momenta
  VectorWaveFunction vr1(rescaledMomenta()[0],
 			 ext[0]->dataPtr(), incoming);
  VectorWaveFunction vr2(rescaledMomenta()[1],
			 ext[1]->dataPtr(), incoming);
  TensorWaveFunction tr3(rescaledMomenta()[2],
			 ext[2]->dataPtr(), outgoing);
  VectorWaveFunction vr4(rescaledMomenta()[3],
			 ext[3]->dataPtr(), outgoing);
  for( unsigned int ihel = 0; ihel < 2; ++ihel ) {
    vr1.reset(2*ihel);
    v1[ihel] = vr1;
    vr2.reset(2*ihel);
    v2[ihel] = vr2;
    tr3.reset(4*ihel);
    t3[4*ihel] = tr3;
    vr4.reset(2*ihel);
    v4[2*ihel] = vr4;
  }
  if( !md ) {
    vr4.reset(1);
    v4[1] = vr4;
  }
  if( !mc ) {
    for(unsigned int ihel=1;ihel<4;++ihel) {
      tr3.reset(ihel);
      t3[ihel] = tr3;
    }
  }
  double dummy(0.);
  ProductionMatrixElement pme = vv2tvHeME(v1, v2, t3, v4,dummy,false);
  createVertex(pme,ext);
}
