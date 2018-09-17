// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoSFFDecayer class.
//

#include "StoSFFDecayer.h"
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

IBPtr StoSFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoSFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoSFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << ten_ << RSfer_ << four_;
}

void StoSFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> ten_ >> RSfer_ >> four_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StoSFFDecayer,GeneralThreeBodyDecayer>
describeHerwigStoSFFDecayer("Herwig::StoSFFDecayer", "Herwig.so");

void StoSFFDecayer::Init() {

  static ClassDocumentation<StoSFFDecayer> documentation
    ("The StoSFFDecayer class implements the general decay of a scalar to "
     "a scalar and two fermions.");

}

WidthCalculatorBasePtr StoSFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<StoSFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		  relativeError()));
}

void StoSFFDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  RSfer_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  four_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    // four point vertex
    if(!offshell) {
      four_[ix] = dynamic_ptr_cast<AbstractFFSSVertexPtr>(current.vertices.first);
      continue;
    }
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractSSSVertexPtr vert1 = dynamic_ptr_cast<AbstractSSSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVSSVertexPtr vert1 = dynamic_ptr_cast<AbstractVSSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractSSTVertexPtr vert1 = dynamic_ptr_cast<AbstractSSTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin3Half) {
      AbstractRFSVertexPtr vert1 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	(current.vertices.first);
      AbstractRFSVertexPtr vert2 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a RS fermion diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      RSfer_[ix] = make_pair(vert1, vert2);
    }
  }
}

double StoSFFDecayer::me2(const int ichan, const Particle & inpart,
			  const ParticleVector & decay,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != inpart.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),
				Helicity::incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true);
    for(unsigned int ix=0;ix<decay.size();++ix) {
      if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
	ScalarWaveFunction::constructSpinInfo(decay[ix],Helicity::outgoing,true);
      }
      else {
	SpinorWaveFunction::
	  constructSpinInfo(outspin_[ix].first,decay[ix],Helicity::outgoing,true);
      }
    }
    return 0.;
  }
  // get the wavefunctions for all the particles
  ScalarWaveFunction outScalar;
  unsigned int isca(0);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
      isca = ix;
      outScalar = ScalarWaveFunction(decay[ix]->momentum(),
				     decay[ix]->dataPtr(),Helicity::outgoing);
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(outspin_[ix].first,decay[ix],Helicity::outgoing);
      outspin_[ix].second.resize(2);
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
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca==0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==2 ? PDT::Spin0 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca == 0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 2 ? PDT::Spin0 : PDT::Spin1Half)));
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  for(unsigned int s1 = 0;s1 < 2; ++s1) {
    for(unsigned int s2 = 0;s2 < 2; ++s2) {
      flows = vector<Complex>(ncf, Complex(0.));
      largeflows = vector<Complex>(ncf, Complex(0.));
      unsigned int idiag(0);
      for(vector<TBDiagram>::const_iterator dit = getProcessInfo().begin();
	  dit != getProcessInfo().end(); ++dit) {
	// channels if selecting
	if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	  ++idiag;
	  continue;
	}
	tcPDPtr offshell = dit->intermediate;
	Complex diag;
	if(offshell) {
	  if(cc&&offshell->CC()) offshell=offshell->CC();
	  double sign = out3[dit->channelType] < out2[dit->channelType] ? 1. : -1.;
	  // intermediate scalar
	  if     (offshell->iSpin() == PDT::Spin0) {
	    ScalarWaveFunction inters = sca_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
	    unsigned int h1(s1),h2(s2);
	    if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	    if(decay[out2[dit->channelType]]->id()<0&&
	       decay[out3[dit->channelType]]->id()>0) {
	      diag =-sign*sca_[idiag].second->
		evaluate(scale,
			 outspin_[out2[dit->channelType]].first [h1],
			 outspin_[out3[dit->channelType]].second[h2],inters);
	    }
	    else {
	      diag = sign*sca_[idiag].second->
		evaluate(scale,
			 outspin_[out3[dit->channelType]].first [h2],
			 outspin_[out2[dit->channelType]].second[h1],inters);
	    }
	  }
	  // intermediate fermion
	  else if(offshell->iSpin() == PDT::Spin1Half) {
	    int iferm = 
	      decay[out2[dit->channelType]]->dataPtr()->iSpin()==PDT::Spin1Half
	      ? out2[dit->channelType] : out3[dit->channelType];
	    unsigned int h1(s1),h2(s2);
	    if(dit->channelType>iferm) swap(h1,h2);
	    sign = iferm<dit->channelType ? 1. : -1.;
	    
	    
	    if((decay[dit->channelType]->id() < 0 &&decay[iferm]->id() > 0 ) ||
	       (decay[dit->channelType]->id()*offshell->id()>0)) {
	      SpinorWaveFunction    inters = fer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].first [h1],swave_);
	      diag = -sign*fer_[idiag].second->
		evaluate(scale,inters,outspin_[iferm].second[h2],outScalar);
	    }
	    else {
	    SpinorBarWaveFunction inters = fer_[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       outspin_[dit->channelType].second[h1],swave_);
	    diag =  sign*fer_[idiag].second->
	      evaluate(scale,outspin_[iferm].first [h2],inters,outScalar);
	    }
	  }
	  // intermediate vector
	  else if(offshell->iSpin() == PDT::Spin1) {
	    VectorWaveFunction interv = vec_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
	    unsigned int h1(s1),h2(s2);
	    if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(decay[out2[dit->channelType]]->id()<0&&
	     decay[out3[dit->channelType]]->id()>0) {
	    diag =-sign*vec_[idiag].second->
	      evaluate(scale,
		       outspin_[out2[dit->channelType]].first [h1],
		       outspin_[out3[dit->channelType]].second[h2],interv);
	  }
	  else {
	    diag = sign*vec_[idiag].second->
	      evaluate(scale,
		       outspin_[out3[dit->channelType]].first [h2],
		       outspin_[out2[dit->channelType]].second[h1],interv);
	  }
	  }
	  // intermediate tensor
	  else if(offshell->iSpin() == PDT::Spin2) {
	    TensorWaveFunction intert = ten_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
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
	  // intermediate RS fermion
	  else if(offshell->iSpin() == PDT::Spin3Half) {
	    int iferm = 
	      decay[out2[dit->channelType]]->dataPtr()->iSpin()==PDT::Spin1Half
	      ? out2[dit->channelType] : out3[dit->channelType];
	    unsigned int h1(s1),h2(s2);
	    if(dit->channelType>iferm) swap(h1,h2);
	    sign = iferm<dit->channelType ? 1. : -1.;
	    if((decay[dit->channelType]->id() < 0 &&decay[iferm]->id() > 0 ) ||
	       (decay[dit->channelType]->id()*offshell->id()>0)) {
	      RSSpinorWaveFunction    inters = RSfer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].first [h1],swave_);
	      diag = -sign*RSfer_[idiag].second->
		evaluate(scale,inters,outspin_[iferm].second[h2],outScalar);
	    }
	    else {
	      RSSpinorBarWaveFunction inters = RSfer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].second[h1],swave_);
	      diag =  sign*RSfer_[idiag].second->
		evaluate(scale,outspin_[iferm].first [h2],inters,outScalar);
	    }
	  }
	  // unknown
	  else throw Exception()
		 << "Unknown intermediate in StoSFFDecayer::me2()" 
		 << Exception::runerror;
	}
	// four point diagram
	else {
	    unsigned int o2 = isca > 0 ? 0 : 1;
	    unsigned int o3 = isca < 2 ? 2 : 1;
	    if(decay[o2]->id() < 0 && decay[o3]->id() > 0) {
	      diag =-four_[idiag]->
	    	evaluate(scale, outspin_[o2].first[s1],
	    		 outspin_[o3].second[s2], outScalar, swave_);
	    }
	    else {
	      diag = four_[idiag]->
	    	evaluate(scale, outspin_[o3].first[s2],
	    		 outspin_[o2].second[s1], outScalar, swave_);
	    }
	}
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
      }
      for(unsigned int ix = 0; ix < ncf; ++ix) {
	if(isca == 0) {
	  (*mes[ix])(0, 0, s1, s2) = flows[ix];
	  (*mel[ix])(0, 0, s1, s2) = largeflows[ix];
	}
	else if(isca == 1 ) { 
	  (*mes[ix])(0, s1, 0, s2) = flows[ix];
	  (*mel[ix])(0, s1, 0, s2) = largeflows[ix];
	}
	else if(isca == 2) { 
	  (*mes[ix])(0, s1,s2, 0) = flows[ix];
	  (*mel[ix])(0, s1,s2, 0) = largeflows[ix] ;
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
