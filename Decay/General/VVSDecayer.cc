// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSDecayer class.
//

#include "VVSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
IBPtr VVSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVSDecayer::fullclone() const {
  return new_ptr(*this);
}
void VVSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    outgoingVertexV_[inter] = AbstractVVVVertexPtr();  
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
      else 
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void VVSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_ << outgoingVertexS_ << outgoingVertexV_;
}

void VVSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_ >> outgoingVertexS_ >> outgoingVertexV_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVSDecayer,GeneralTwoBodyDecayer>
describeHerwigVVSDecayer("Herwig::VVSDecayer", "Herwig.so");

void VVSDecayer::Init() {

  static ClassDocumentation<VVSDecayer> documentation
    ("The VVSDecayer class implements the decay of a vector"
     " to a vector and a scalar");

}

double VVSDecayer::me2(const int , const Particle & inpart,
 		       const ParticleVector & decay,
		       MEOption meopt) const {
  bool massless = ( decay[0]->id()==ParticleID::gamma || 
		    decay[0]->id()==ParticleID::g );
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors_[1],decay[0],outgoing,true,massless);
    ScalarWaveFunction::
      constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::
    calculateWaveFunctions(vectors_[1],decay[0],outgoing,massless);
  ScalarWaveFunction sca(decay[1]->momentum(),decay[1]->dataPtr(),outgoing);
  Energy2 scale(sqr(inpart.mass()));
  for(unsigned int in=0;in<3;++in) {
    for(unsigned int out=0;out<3;++out) {
      if(massless&&out==1) ++out;
      (*ME())(in,out,0) = 0.;
      for(auto vert : vertex_)
	(*ME())(in,out,0) += 
	  vert->evaluate(scale,vectors_[0][in],vectors_[1][out],sca);
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(inpart.dataPtr(),decay[0]->dataPtr(),
			 decay[1]->dataPtr());
  // return the answer
  return output;
}

Energy VVSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outb.first->iSpin() == PDT::Spin0 )
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, 
				       outa.first, outb.first);
    else {
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, 
				       outb.first, outa.first);
      swap(mu1sq, mu2sq);
    }
    double vn = norm(perturbativeVertex_[0]->norm());
    if(vn == ZERO || mu1sq == ZERO) return ZERO;
    double me2 = 2. + 0.25*sqr(1. + mu1sq - mu2sq)/mu1sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = vn*me2*pcm/(24.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VVSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  unsigned int ivec(0),isca(1);
  if(decay[ivec]->dataPtr()->iSpin()!=PDT::Spin1) swap(ivec,isca);
  if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vectors3_[0], rho3_,
					       const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vectors3_[0] ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[1],decay[ivec],outgoing,true,false); 
    ScalarWaveFunction::
      constructSpinInfo(decay[isca],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_      ,decay[2],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin1,
								       PDT::Spin0, PDT::Spin1)));
  bool massless= decay[ivec]->mass()!=ZERO;
  // create wavefunctions
  VectorWaveFunction::calculateWaveFunctions(vectors3_[1],decay[0],outgoing,massless);
  ScalarWaveFunction scal(decay[isca]->momentum(),  decay[isca]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_      ,decay[2],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[2]->momentum(),
  					  decay[2]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  Energy2 scale(sqr(inpart.mass()));
  
  const GeneralTwoBodyDecayer::CFlow & colourFlow
    = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif

   for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
     for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
       if(massless && iv1==1) continue;
       for(unsigned int ig = 0; ig < 2; ++ig) {
	 // radiation from the incoming vector
	 if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(incomingVertex_[inter]);
	   VectorWaveFunction vectorInter = 
	     incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vectors3_[0][iv0],
					      gluon_[2*ig],inpart.mass());
	   
	   assert(vectors3_[0][iv0].particle()->id()==vectorInter.particle()->id());
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectorInter,vectors3_[1][iv1],scal);
	   if(!couplingSet) {
	     gs = abs(incomingVertex_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	     (*ME[colourFlow[0][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[0][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the outgoing vector
	 if((decay[ivec]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[ivec]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertexV_[inter]);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[ivec]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],vectors3_[1][iv1],decay[ivec]->mass());
	   
	   assert(vectors3_[1][iv1].particle()->id()==vectorInter.particle()->id());
	   
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv0],vectorInter,scal);
	   if(!couplingSet) {
	     gs = abs(outgoingVertexV_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	     (*ME[colourFlow[1][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[1][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the outgoing scalar
	 if((decay[isca]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[isca]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertexS_[inter]);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[isca]->dataPtr();
	   if(off->CC()) off = off->CC();
	   ScalarWaveFunction scalarInter = 
	     outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],scal,decay[isca]->mass());
	
	   assert(scal.particle()->id()==scalarInter.particle()->id());
	 
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv0],vectors3_[1][iv1],scalarInter);
	   if(!couplingSet) {
	     gs = abs(outgoingVertexS_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	     (*ME[colourFlow[2][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[2][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
       }
     }
   }
   
   // contract matrices 
   double output=0.;
   for(unsigned int ix=0; ix<nflow; ++ix){
     for(unsigned int iy=0; iy<nflow; ++iy){
       output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
     }
   }
   // divide by alpha_(S,EM)
   output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
   double ratio = output/total;
   if(abs(ratio)>1e-20) {
     generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
     for(unsigned int ix=0;ix<decay.size();++ix)
       generator()->log() << *decay[ix] << "\n";
     generator()->log() << "Test of gauge invariance " << ratio << "\n";
   }
#endif
   // return the answer
   return output;
}
