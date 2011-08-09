// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoFFFFDecayer class.
//

#include "StoFFFFDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include <numeric>

using namespace Herwig;

IBPtr StoFFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoFFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoFFFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << firstVVS_ << firstVSS_ << firstSSS_ << firstFFS_
     << secondFFV_ << secondFFS_
     << thirdFFV_ << thirdFFS_ << sign_;
}

void StoFFFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> firstVVS_ >> firstVSS_ >> firstSSS_ >> firstFFS_
     >> secondFFV_ >> secondFFS_
     >> thirdFFV_ >> thirdFFS_ >> sign_;
}

DescribeClass<StoFFFFDecayer,GeneralFourBodyDecayer>
describeStoFFFFDecayer("Herwig::StoFFFFDecayer", "Herwig.so");

void StoFFFFDecayer::Init() {

  static ClassDocumentation<StoFFFFDecayer> documentation
    ("The StoFFFFDecayer class performs the 4-fermion "
     "decays of scalar particles in BSM models");

}

void StoFFFFDecayer::doinit() {
  GeneralFourBodyDecayer::doinit(); 
  unsigned int ndiags = getProcessInfo().size();
  firstVVS_ .resize(ndiags);
  firstVSS_ .resize(ndiags);
  firstSSS_ .resize(ndiags);
  firstFFS_ .resize(ndiags);
  secondFFV_.resize(ndiags);
  secondFFS_.resize(ndiags);  
  thirdFFV_ .resize(ndiags);
  thirdFFS_ .resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    const NBDiagram & current = getProcessInfo()[ix];
    // first vertex
    firstVVS_[ix] = dynamic_ptr_cast<AbstractVVSVertexPtr>(current.vertex);
    firstVSS_[ix] = dynamic_ptr_cast<AbstractVSSVertexPtr>(current.vertex);
    firstSSS_[ix] = dynamic_ptr_cast<AbstractSSSVertexPtr>(current.vertex);
    firstFFS_[ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(current.vertex);
    // get the other vertices
    const NBVertex & second = current.vertices.begin()->second.incoming ?
      current.vertices.begin()->second : (++current.vertices.begin())->second;
    // get the other vertices
    const NBVertex & third = current.vertices.begin()->second.incoming ?
      (++current.vertices.begin())->second : (++second.vertices.begin())->second;
    // second vertex
    secondFFV_[ix] = dynamic_ptr_cast<AbstractFFVVertexPtr>(second.vertex);
    secondFFS_[ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(second.vertex);  
    // third vertex
    thirdFFV_ [ix] = dynamic_ptr_cast<AbstractFFVVertexPtr>(third .vertex);
    thirdFFS_ [ix] = dynamic_ptr_cast<AbstractFFSVertexPtr>(third .vertex); 
    assert( ( firstVVS_[ix] ||  firstVSS_[ix] || 
	      firstSSS_[ix] ||  firstFFS_[ix] ) && 
	    (secondFFV_[ix] || secondFFS_[ix] ) &&
	    ( thirdFFV_[ix] ||  thirdFFS_[ix] ));
    // NO sign
    switch(current.channelType) {
    case 1234: case 1342: case 1423: 
    case 2143: case 2314: case 2431:
    case 3124: case 3241: case 3412:
    case 4132: case 4213: case 4321:
      sign_.push_back( 1.);
      break;
    case 1243: case 1324: case 1432:
    case 2134: case 2341: case 2413:
    case 3142: case 3214: case 3421:
    case 4123: case 4231: case 4312:
      sign_.push_back(-1.);
      break;
    default:
      assert(false);
    }
  }
}

double StoFFFFDecayer::me2(const int ichan, const Particle & inpart,
			   const ParticleVector & decay, MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming->id() != inpart.id();
  // special handling or first/last call
  const vector<vector<double> > & cfactors(getColourFactors());
  const vector<vector<double> > & nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  Energy2 scale(sqr(inpart.mass()));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&inpart),
			     Helicity::incoming);
    swave_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),
				Helicity::incoming);
  }
  // setup spin info when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
			Helicity::incoming,true);
    // outgoing particles
    for(unsigned int ix = 0; ix < 4; ++ix) {
      SpinorWaveFunction::
      constructSpinInfo(outwave_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 4; ++ix) {
    SpinorWaveFunction::
      calculateWaveFunctions(outwave_[ix].first,decay[ix],Helicity::outgoing);
    outwave_[ix].second.resize(2);
    if(outwave_[ix].first[0].wave().Type() == u_spinortype) {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].first[iy].conjugate();
      }
    }
    else {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].second[iy].conjugate();
      }
    }
  }
  // matrix element for the colour flows
  vector<Complex> flows(ncf, Complex(0.)),largeflows(ncf, Complex(0.));
  vector<DecayMatrixElement> mes(ncf,DecayMatrixElement(PDT::Spin0,
							PDT::Spin1Half,PDT::Spin1Half,
							PDT::Spin1Half,PDT::Spin1Half));
  vector<DecayMatrixElement> mel(ncf,DecayMatrixElement(PDT::Spin0,
							PDT::Spin1Half,PDT::Spin1Half,
							PDT::Spin1Half,PDT::Spin1Half));
  // calculate the matrix element
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag=0;
	  for(vector<NBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    if(ichan>=0&&diagramMap()[ichan]!=idiag) {
	      ++idiag;
	      continue;
	    }
	    // location of the particles
	    int iloc[4]={ dit->channelType/1000     -1,
			 (dit->channelType%1000)/100-1,
			 (dit->channelType%100)/10  -1,
			  dit->channelType%10       -1};
	    // NO sign
	    double sign = sign_[idiag];
	    // work out the type of topology
	    bool topo = dit->vertices.begin()->second.incoming;
	    const NBVertex & second = topo ? 
   	          dit->vertices.begin()  ->second :
	      (++(dit->vertices.begin()))->second;
	    const NBVertex & third = topo ?
	      (++(dit->  vertices.begin()))->second :
	      (++(second.vertices.begin()))->second;
	    // extract the intermediates
	    tPDPair inter = make_pair(second.incoming,
				      third .incoming->CC());
	    if(cc&&inter.first ->CC()) inter.first  = inter.first ->CC();
	    if(cc&&inter.second->CC()) inter.second = inter.second->CC();

// 	    //\todo remove testing
// 	    if(abs(inter.first ->id())!=ParticleID::t||
// 	       abs(inter.second->id())!=ParticleID::Wplus) {
// 	      ++idiag;
// 	      continue;
// 	    }

	    // value for the diagram
	    Complex diag(0.);
	    // first compute the last part of the diagram
	    VectorWaveFunction offVector2;
	    ScalarWaveFunction offScalar2;
	    // intermeidate scalar
	    if(inter.second->iSpin()==PDT::Spin0) {
	      if(decay[iloc[2]]->id()<0&&
		 decay[iloc[3]]->id()>0) {
		sign *= -1.;
		offScalar2 = thirdFFS_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[2]].first [ihel[iloc[2]]],
			   outwave_[iloc[3]].second[ihel[iloc[3]]]);
	      }
	      else {
		offScalar2 = thirdFFS_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[3]].first [ihel[iloc[3]]],
			   outwave_[iloc[2]].second[ihel[iloc[2]]]);
	      }
	    }
	    // intermediate vector
	    else if(inter.second->iSpin()==PDT::Spin1) {
	      if(decay[iloc[2]]->id()<0&&
		 decay[iloc[3]]->id()>0) {
		sign *= -1.;
		offVector2 = thirdFFV_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[2]].first [ihel[iloc[2]]],
			   outwave_[iloc[3]].second[ihel[iloc[3]]]);
	      }
	      else {
		offVector2 = thirdFFV_[idiag]->
		  evaluate(scale, widthOption(),inter.second,
			   outwave_[iloc[3]].first [ihel[iloc[3]]],
			   outwave_[iloc[2]].second[ihel[iloc[2]]]);
	      }
	    }
	    // unknown
	    else
	      assert(false);
	    // first topology
	    if(topo) {
	      // first intermediate
	      if(inter.first->CC()) inter.first = inter.first->CC();
	      VectorWaveFunction offVector1;
	      ScalarWaveFunction offScalar1;
	      // intermeidate scalar
	      if(inter.first->iSpin()==PDT::Spin0) {
		if(decay[iloc[0]]->id()<0&&
		   decay[iloc[1]]->id()>0) {
		  sign *= -1.;
		  offScalar1 = secondFFS_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[0]].first [ihel[iloc[0]]],
			     outwave_[iloc[1]].second[ihel[iloc[1]]]);
		}
		else {
		  offScalar1 = secondFFS_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[1]].first [ihel[iloc[1]]],
			     outwave_[iloc[0]].second[ihel[iloc[0]]]);
		}
	      }
	      // intermediate vector
	      else if(inter.first->iSpin()==PDT::Spin1) {
		if(decay[iloc[0]]->id()<0&&
		   decay[iloc[1]]->id()>0) {
		  sign *= -1.;
		  offVector1 = secondFFV_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[0]].first [ihel[iloc[0]]],
			     outwave_[iloc[1]].second[ihel[iloc[1]]]);
		}
		else {
		  offVector1 = secondFFV_[idiag]->
		    evaluate(scale, widthOption(),inter.first,
			     outwave_[iloc[1]].first [ihel[iloc[1]]],
			     outwave_[iloc[0]].second[ihel[iloc[0]]]);
		}
	      }
	      // unknown
	      else
		assert(false);
	      // put it all together	      
	      if(inter.first->iSpin()==PDT::Spin0) {	      
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = firstSSS_[idiag]->
		    evaluate(scale,swave_,offScalar1,offScalar2);
		}
		else if(inter.second->iSpin()==PDT::Spin1) {
		  diag = firstVSS_[idiag]->
		    evaluate(scale,offVector2,offScalar1,swave_);
		}
	      }
	      else if(inter.first->iSpin()==PDT::Spin1) {	      
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = firstVSS_[idiag]->
		    evaluate(scale,offVector1,offScalar2,swave_);
		}
		else if(inter.second->iSpin()==PDT::Spin1) {
		  diag = firstVVS_[idiag]->
		    evaluate(scale,offVector1,offVector2,swave_);
		}
	      }
	    }
	    // second topology
	    else {
	      if(decay[iloc[0]]->id()<0&&decay[iloc[1]]->id()>0) {
		sign *= -1.;
		SpinorWaveFunction    inters = firstFFS_[idiag]->
		  evaluate(scale,widthOption(),inter.first,
			   outwave_[iloc[0]].first [ihel[iloc[0]]],swave_);
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = secondFFS_[idiag]->
		    evaluate(scale,inters,outwave_[iloc[1]].second[ihel[iloc[1]]],
			     offScalar2);
		}
		else {
		  diag = secondFFV_[idiag]->
		    evaluate(scale,inters,outwave_[iloc[1]].second[ihel[iloc[1]]],
			     offVector2);
		}
	      }
	      else {
		SpinorBarWaveFunction inters = firstFFS_[idiag]->
		  evaluate(scale,widthOption(),inter.first,
			   outwave_[iloc[0]].second[ihel[iloc[0]]],swave_);
		if(inter.second->iSpin()==PDT::Spin0) {
		  diag = secondFFS_[idiag]->
		    evaluate(scale,outwave_[iloc[1]].first [ihel[iloc[1]]],inters,
			     offScalar2);
		}
		else {
		  diag = secondFFV_[idiag]->
		    evaluate(scale,outwave_[iloc[1]].first [ihel[iloc[1]]],inters,
			     offVector2);
		}
	      }
	    }
 	    // apply NO sign
 	    diag *= sign;
	    // matrix element for the different colour flows
	    if(ichan<0) {
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
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
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
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    mes[ix](0,ihel[0],ihel[1],ihel[2],ihel[3]) =      flows[ix];
	    mel[ix](0,ihel[0],ihel[1],ihel[2],ihel[3]) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix].contract(mes[iy],rho_)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix].contract(mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal-=pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow].contract(mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2*scale*UnitRemoval::InvE2;
}
