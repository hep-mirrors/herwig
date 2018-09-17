// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VtoFFVDecayer class.
//

#include "VtoFFVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr VtoFFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VtoFFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void VtoFFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << ten_;
}

void VtoFFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> ten_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VtoFFVDecayer,GeneralThreeBodyDecayer>
describeHerwigVtoFFVDecayer("Herwig::VtoFFVDecayer", "Herwig.so");

void VtoFFVDecayer::Init() {

  static ClassDocumentation<VtoFFVDecayer> documentation
    ("The VtoFFVDecayer class implements the general three-body "
     "decay of a vector to a two fermions and a vector.");

}

WidthCalculatorBasePtr VtoFFVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<VtoFFVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass(),relativeError()));
}

void VtoFFVDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVVVertexPtr vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractVVTVertexPtr vert1 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
  }
}

double VtoFFVDecayer::me2(const int ichan, const Particle & inpart,
			  const ParticleVector & decay,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    VectorWaveFunction::
      calculateWaveFunctions(inVector_,rho_,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }

  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(inVector_,const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true,false);
    for(unsigned int ix=0;ix<decay.size();++ix) {
      if(decay[ix]->dataPtr()->iSpin()==PDT::Spin1) {
	VectorWaveFunction::constructSpinInfo(outVector_,decay[ix],
					      Helicity::outgoing,true,false);
      }
      else {
	SpinorWaveFunction::
	  constructSpinInfo(outspin_[ix].first,decay[ix],Helicity::outgoing,true);
      }
    }
  }
  unsigned int ivec(0);
  bool massless = false;
  for(unsigned int ix = 0; ix < decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin() == PDT::Spin1) {
      massless = decay[ix]->id()==ParticleID::g || decay[ix]->id()==ParticleID::gamma;
      ivec = ix;
      VectorWaveFunction::
	calculateWaveFunctions(outVector_, decay[ix], Helicity::outgoing, massless );
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(outspin_[ix].first,decay[ix],Helicity::outgoing);
      outspin_[ix].second.resize(2);
      // Need a ubar and a v spinor
      if(outspin_[ix].first[0].wave().Type() == SpinorType::u) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(inpart.mass()));
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.));
  // setup the DecayMatrixElement
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  
  //the channel possiblities
  static const unsigned int out2[3] = {1,0,0}, out3[3] = {2,2,1};
  for(unsigned int vi = 0; vi < 3; ++vi) {
    for(unsigned int s1 = 0; s1 < 2; ++s1) {
      for(unsigned int s2 = 0; s2 < 2; ++s2) {
	for(unsigned int v1 = 0; v1 < 3; ++v1) {
	  if ( massless && v1 == 1 ) continue; 
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag(0);
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    // channels if selecting
	    if( ichan >=0 && diagramMap()[ichan] != idiag ) {
	      ++idiag;
	      continue;
	    }
	    tcPDPtr offshell = dit->intermediate;
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    Complex diag;
	    unsigned int o2(out2[dit->channelType]), o3(out3[dit->channelType]);
	    double sign = (o3 < o2) ? 1. : -1.;
	    // intermediate scalar
	    if(offshell->iSpin() == PDT::Spin0) {
	      ScalarWaveFunction inters = sca_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], 
			 inVector_[vi]);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1, h2);
	      if(decay[o2]->id() < 0 &&  decay[o3]->id() > 0) {
		diag = -sign*sca_[idiag].second->
		  evaluate(scale,outspin_[o2].first[h1],
			   outspin_[o3].second[h2],inters);
	      }
	      else {
		diag = sign*sca_[idiag].second->
		  evaluate(scale, outspin_[o3].first [h2],
			   outspin_[o2].second[h1],inters);
	      }
	    }
	    // intermediate fermion
	    else if(offshell->iSpin() == PDT::Spin1Half) {
	      int iferm = (decay[o2]->dataPtr()->iSpin() == PDT::Spin1Half) 
		? o2 : o3;
	      unsigned int h1(s1),h2(s2);
	      if(dit->channelType > iferm) swap(h1, h2);
	      sign = iferm < dit->channelType ? 1. : -1.;
	      if(decay[dit->channelType]->id() < 0 && decay[iferm]->id() > 0 ) {
		SpinorWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].first[h1], inVector_[vi]);
		diag = -sign*fer_[idiag].second->
		  evaluate(scale,inters,outspin_[iferm].second[h2], outVector_[v1]);
	      }
	      else {
		SpinorBarWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].second[h1],inVector_[vi]);
		diag =  sign*fer_[idiag].second->
		  evaluate(scale,outspin_[iferm].first [h2],inters, outVector_[v1]);
	      }
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = vec_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], 
			 inVector_[vi]);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1,h2);
	      if(decay[o2]->id() < 0 && decay[o3]->id() > 0) {
		diag =-sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o2].first[h1],
			   outspin_[o3].second[h2], interv);
	      }
	      else {
		diag = sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o3].first[h2],
			   outspin_[o2].second[h1], interv);
	      }
	    }
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = ten_[idiag].first->
		evaluate(scale, widthOption(), offshell, inVector_[vi], 
			 outVector_[v1]);
	      unsigned int h1(s1),h2(s2);
	      if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	      if(decay[out2[dit->channelType]]->id()<0&&
		 decay[out3[dit->channelType]]->id()>0) {
		diag =-sign*ten_[idiag].second->
		  evaluate(scale,
			   outspin_[out2[dit->channelType]].first [h1],
			   outspin_[out3[dit->channelType]].second[h2],intert);
	      }
	      else {
		diag = sign*ten_[idiag].second->
		  evaluate(scale,
			   outspin_[out3[dit->channelType]].first [h2],
			   outspin_[out2[dit->channelType]].second[h1],intert);
	      }
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in VtoFFVDecayer::me2()" 
	      << Exception::runerror;
	    
	    // matrix element for the different colour flows
	    if(ichan < 0) {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    else {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }	    
	    ++idiag;
	  } //end of diagrams
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    if (ivec == 0) {
	      (*mes[ix])(vi, v1, s1, s2) = flows[ix];
	      (*mel[ix])(vi, v1, s1, s2) = largeflows[ix];
	    }
	    else if(ivec == 1) {
	      (*mes[ix])(vi, s1, v1, s2) = flows[ix];
	      (*mel[ix])(vi, s1, v1, s2) = largeflows[ix];
	      
	    }
	    else if(ivec == 2) { 
	      (*mes[ix])(vi, s1, s2, v1) = flows[ix];
	      (*mel[ix])(vi, s1, s2, v1) = largeflows[ix];
	    }
	  }

	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}
