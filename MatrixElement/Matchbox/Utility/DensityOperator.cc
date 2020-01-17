// -*- C++ -*-
//
// DensityOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DensityOperator class.
//

#include "DensityOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXCombData.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;


DensityOperator::DensityOperator() : Nc(3.0), TR(0.5) { }

DensityOperator::~DensityOperator() {}

void DensityOperator::clear() {
	theCorrelatorMap.clear();
}

// Prepare density operator if not done before
void DensityOperator::prepare(const cPDVector& mePartonData) {

	// Take parton data and create key of type 33bar888 to use as key for basis
	vector<PDT::Colour> mePartonDataColoured = theColourBasis->normalOrderMap(mePartonData);

	if( theDensityOperatorMap.count( mePartonDataColoured ) == 0 ){
		// Get basis dimension
		size_t dim = theColourBasis->prepare( mePartonData, false );
		// Allocate space for density matrix
		theDensityOperatorMap.insert(make_pair(mePartonDataColoured,matrix<Complex> (dim,dim)));
	}
}

// Fill the density matrix for the first time
void DensityOperator::fill(const Ptr<MatchboxXComb>::ptr MBXCombPtr,
				 const cPDVector& partons,
				 const vector<Lorentz5Momentum>& momenta){

	// Map from helicity structure to amplitude vector in the color basis
	const map<vector<int>,CVector>& amplitudeMap = MBXCombPtr->lastAmplitudes();
	const cPDVector& mePartonData = MBXCombPtr->mePartonData();
	// Get the dimension of the basis for the indexChange method
	size_t dim = theColourBasis->prepare(mePartonData,false);
	// Normal order partons according to ColourBasis
	const vector<PDT::Colour> mePartonDataColoured = theColourBasis->normalOrderMap(mePartonData);

	// Get the colour basis to colour basis map for this hard subprocess
	map<cPDVector,map<size_t,size_t> >::iterator cb2cbit = theColourBasisToColourBasisMap.find(mePartonData);
	// Fill the map with this hard subprocess if it doesn't have it yet
	if ( cb2cbit == theColourBasisToColourBasisMap.end() ) {
		// Get the index map for the partons as they are ordered in the MatchboxXComb
		// object
		const map<size_t,size_t>& mbIndexMap = theColourBasis->indexMap().at(mePartonData);
		
		// Get the index map for the partons as they are ordered
		// in the shower.
		// Use prepare, as it might not have been done for this order of the partons
		theColourBasis->prepare(partons,false);
		const map<size_t,size_t>& showerIndexMap = theColourBasis->indexMap().at(partons);
		
		// ensure the maps are of the same size
		assert( mbIndexMap.size() == showerIndexMap.size() );

		// Loop over both sets of partons to determine how the order between them
		// differs, then translate the index to the colour basis and put the
		// key-value pair into the map
		map<size_t,size_t> cb2cbMap;

		// Get the momenta for comparison
		const vector<Lorentz5Momentum>& meMomenta = MBXCombPtr->matchboxME()->lastMEMomenta();
		// Make sure there's the same number of momenta in both vectors
		assert( momenta.size() == meMomenta.size() );
		bool done;
		// Boost the momenta to the same frame
		const vector<Lorentz5Momentum> momentaCM = boostToRestFrame(momenta);
		for ( size_t i = 0; i < momentaCM.size(); i++ ) {
		  // The cb2cb map is intended to translate how the indices
		  // are different in the colour basis due to the different 
		  // orderings of the partons, so it should only bother
		  // with coloured particles.
                 if ( partons[i]->coloured() ) {
			done = false;
			for ( size_t j = 0; j < meMomenta.size(); j++ ) {
	if ( !done ) {
		if ( compareMomentum(momentaCM[i],meMomenta[j]) ) {
				cb2cbMap[mbIndexMap.at(i)] = showerIndexMap.at(j);
				done = true;
		}
	}	
			}
		 }
		}

		// Make sure all momenta have been identified
		assert( cb2cbMap.size() == mePartonDataColoured.size() );

		// Add the map to the cache
		theColourBasisToColourBasisMap[mePartonData] = cb2cbMap;

		// Get the iterator
		cb2cbit = theColourBasisToColourBasisMap.find(mePartonData);
	}

	// With the cb2cb index map we can create the basis vector index map
	const map<size_t,size_t> vectorMap = theColourBasis->indexChange(mePartonDataColoured,
									 dim,cb2cbit->second);
	
	// Prepare density operator (allocate) for set of particles
	prepare(mePartonData);

	// Check that density operator (place holder for) exist in density operator map
	map<vector<PDT::Colour>,matrix<Complex> >::iterator dOit = theDensityOperatorMap.find(mePartonDataColoured);
	assert(dOit != theDensityOperatorMap.end());

	// Initialize the density operator
	matrix<Complex>& densOp = dOit->second;
	for(unsigned int i = 0; i < (amplitudeMap.begin()->second).size(); i++){
		for(unsigned int j = 0; j < (amplitudeMap.begin()->second).size(); j++){
			densOp (i,j) = 0;
		}
	}

	// Fill the density operator, ok to sum over helicities since the density operator
	// exist at the probability level
	CVector amplitude;
	for ( map<vector<int>,CVector>::const_iterator itHel = amplitudeMap.begin();
	itHel != amplitudeMap.end(); itHel++ ) {
		amplitude = itHel->second;
		for ( unsigned int i = 0; i < amplitude.size(); i++ ) {
			for ( unsigned int j = 0; j < amplitude.size(); j++ ) {
	// vectorMap is used such that densOp is filled according to the
	// basis order defined by the input cPDVector partons.
	densOp (vectorMap.at(i),vectorMap.at(j)) += amplitude(i)*std::conj(amplitude(j));
			}
		}
	}
	
	// Colour conservation check
	colourConservation(mePartonData);
}

// Update density matrix after emitting or splitting a gluon
// The emissionsMap argument contains the relation of the indices in the smaller (before)
// and larger basis (after). the first 3-tuple contains the indices
// (emitter before, emitter after, emitted parton)
// The map contains the old and new indices of all other partons (not involved)
void DensityOperator::evolve(const map<pair<size_t,size_t>,Complex>& Vijk, 
					 const cPDVector& before,
					 const cPDVector& after,
					 const map<std::tuple<size_t,size_t,size_t>,map<size_t,size_t> >& emissionsMap,
					 const bool splitAGluon,
					 const bool initialGluonSplitting) {

	size_t dimBefore = theColourBasis->prepare( before, false );
	size_t dimAfter = theColourBasis->prepare( after, false );

	vector<PDT::Colour> beforeColoured = theColourBasis->normalOrderMap(before);
	vector<PDT::Colour> afterColoured = theColourBasis->normalOrderMap(after);

	const map<vector<PDT::Colour>,matrix<Complex> >::iterator dOit = theDensityOperatorMap.find(beforeColoured);
	assert(dOit != theDensityOperatorMap.end());

	const matrix<Complex>& densOpBefore = dOit->second;

	prepare(after);
	matrix<Complex>& densOpAfter = theDensityOperatorMap[afterColoured];
	for(size_t i = 0; i < densOpAfter.size1(); i++){
		for(size_t j = 0; j < densOpAfter.size2(); j++){
			densOpAfter(i,j) = 0;
		}
	}

	compressed_matrix<double> Tij;
	matrix<Complex> TijMn (dimAfter,dimBefore);
	compressed_matrix<double> Tk;
	matrix<Complex> TijMnTkdagger (dimAfter,dimAfter);
	Complex V;
	// Compensate for sign from updated density matrix
	// TODO Check signs again
	double sign = -1.0;
	if ( splitAGluon )
		sign = 1.0;
	// Loop over emitter legs ij and recoil legs k and add the contribution,
	// for Vijk is assumed to contain the factor 4*pi*\alpha_s/pi.pj
	typedef map<std::tuple<size_t,size_t,size_t>,map<size_t,size_t> > dictMap;
	int ij,k;
	
	// Loop over emitters
	for(dictMap::const_iterator ijit = emissionsMap.begin();
			ijit != emissionsMap.end(); ijit++) {
		// get first element in 3-tuple, i.e., emitter before
		ij = std::get<0>(ijit->first);
		assert(before[ij]->coloured());
		// Get rectangular matrices T_{ij} or S_{ij} taking us from the
		// smaller to the larger basis depending on
		// the involved particles before and after, the emitter index before ij,
		// the emitter after and the emitted parton index
		int i=std::get<1>(ijit->first); // emitter after
		int j=std::get<2>(ijit->first);// emitted parton
		Tij = theColourBasis->charge(before,after,
				 ij,i,j,ijit->second).first;
		TijMn = prodSparseDense(Tij,densOpBefore);

		// Loop over spectators
		for(dictMap::const_iterator kit = emissionsMap.begin();
	kit != emissionsMap.end(); kit++) {
			k = std::get<0>(kit->first);
			assert(before[k]->coloured());
			// Standard case of gluon radiation
			if ( ijit != kit || splitAGluon || initialGluonSplitting ) {
	int k_after=std::get<1>(kit->first); // For color structure k now has role of emitter
	int k_emission= std::get<2>(kit->first); // Emitted parton index
	Tk = theColourBasis->charge(before,after,
						k, k_after, k_emission, kit->second).first;
	TijMnTkdagger = prodDenseSparse(TijMn,Tk);
	// sign == -1.0 if it isn't a gluon splitting into a qqbar pair
	V = sign*(1.0/colourNorm(before[ij]))*
		Vijk.at(make_pair(ij,k));
	densOpAfter += V*TijMnTkdagger;
			}
		}
	}
	// Check that the density operator does not vanish
	assert( theColourBasis->me2(after,densOpAfter) != 0.0 );
	colourConservation(after);
}

// The 3-tuple contains (emitter index, spectator index, emission pid)
double DensityOperator::colourMatrixElementCorrection(const std::tuple<size_t,size_t,long>& ikemission,
									const cPDVector& particles) {
	const int i = std::get<0>(ikemission);
	const int k = std::get<1>(ikemission);
	//const long emissionID = std::get<2>(ikemission);


	// Get the density operator
	// normal order particles as in ColourBasis
	vector<PDT::Colour> particlesColoured = theColourBasis->normalOrderMap(particles);
	// ... and find corresponding density operator
	const map<vector<PDT::Colour>,matrix<Complex> >::iterator particlesit = theDensityOperatorMap.find(particlesColoured);
	assert(particlesit != theDensityOperatorMap.end());
	const matrix<Complex>& densOp = particlesit->second;

	double Ti2 = colourNorm(particles[i]);

	// Emitter-spectator pair
	const pair<size_t,size_t> ik = make_pair(i,k);
	// Create key for color matrix element correction map
	pair<vector<PDT::Colour>,pair<size_t,size_t> > particlesAndLegs = make_pair(particlesColoured,ik);

	// Result
	double res = 0;

	// Check if it has already been calculated
	// TODO: move this check earlier (we only need particlesColoured to check if it has been
	//			 calculated).
	// Check for color matrix element associated with key, and calculate if not done
	const map<pair<vector<PDT::Colour>,pair<size_t,size_t> >,double >::const_iterator corrit = theCorrelatorMap.find(particlesAndLegs);
	if ( corrit == theCorrelatorMap.end() ) {
		double corrME2 = theColourBasis->colourCorrelatedME2(ik,particles,densOp);
		double me2 = theColourBasis->me2(particles,densOp);
		res = -(1/Ti2)*corrME2/me2;

		if ( particles[i]->id() == ParticleID::g )
			res *= 2.;
		
		theCorrelatorMap.insert(make_pair(particlesAndLegs,res));

	} else {
		res = corrit->second;
	}

	return res;
}

double DensityOperator::colourNorm(const cPDPtr particle) {
	if ( particle->id() == ParticleID::g ) {
		return Nc; //is 3.0 for Nc = 3, TR = 1/2
	} else if ( particle->iColour() == PDT::Colour3 || particle->iColour() == PDT::Colour3bar ) {
	  return TR*(Nc*Nc-1.)/Nc; // is 4.0/3.0 for Nc = 3, TR = 1/2
	} else {
		throw Exception() << "Colour matrix element corrections only work "
					<< "on quark and gluon legs.	"
					<< Exception::runerror;
	}
}

void DensityOperator::colourConservation(const cPDVector& particles) {
	// To contain (emitter, spectator, emission pid)
	std::tuple<size_t,size_t,long> ikemission;
	// Normal order particles as defined in ColourBasis
	const vector<PDT::Colour> particlesColoured = theColourBasis->normalOrderMap(particles);
	vector<double> sum(particlesColoured.size(),0.0);
	size_t iterm = 0;

	// To compensate for the CMEC having a 1/(1+\delta(i is gluon)) factor
	// in the splitting kernel
	double gluonFactor = 1.0;
	// Loop over "emitters" to check color conservation for
	for ( size_t i = 0; i < particles.size(); i++ ) {
		if ( particles[i]->coloured() ) {
			for ( size_t k = 0; k < particles.size(); k++ ) {
	if ( particles[k]->coloured() && i != k ) {
		ikemission = std::make_tuple(i,k,ParticleID::g);
		if ( particles[i]->id() != ParticleID::g ) {
			gluonFactor = 1.0;
		} else {
			gluonFactor = 1./2.;
		}
		sum[iterm] += gluonFactor*colourMatrixElementCorrection(ikemission,particles);
	}
			}
			iterm++;
		}
	}
	for ( size_t i = 0; i < sum.size(); i++ )
		assert( std::abs(sum[i]-1.0) < pow(10.0,-10.0));
}

matrix<Complex> DensityOperator::prodSparseDense(const compressed_matrix<double>& Tij,
						 const matrix<Complex>& Mn){
	// Dimension before emission
	size_t dimBefore = Tij.size2();
	// Dimension after emission
	size_t dimAfter = Tij.size1();

	//Check matrix dimensions
	assert( dimBefore == Mn.size1() );

	// Allocate memory for the matrix
	matrix<Complex> TijMn (dimAfter,dimBefore,0);
	// Use iterators for the compressed matrix to iterate only over the non-zero
	// elements.
	size_t ii;
	size_t jj;
	for ( compressed_matrix<double>::const_iterator1 it1 = Tij.begin1();
	it1 != Tij.end1(); it1++ ) {
		for ( compressed_matrix<double>::const_iterator2 it2 = it1.begin();
		it2 != it1.end(); it2++ ) {
			ii = it2.index1();
			jj = it2.index2();
			for ( size_t kk = 0; kk < dimBefore; kk++ ) {
	// *it2 is Tij(ii,jj)
	TijMn(ii,kk) += (*it2)*Mn(jj,kk);
			}
		}
	}
	return TijMn;
}

matrix<Complex> DensityOperator::prodDenseSparse(const matrix<Complex>& TijMn,
						 const compressed_matrix<double>& Tk){
	// The compressed matrix comes from the charge method, do not transpose yet
	// Dimension after emission
	size_t dimAfter = Tk.size1();//Since this method returns TijMn*Tk^\dagger

	//Check matrix dimensions
	assert( TijMn.size2() == Tk.size2() );

	// Allocate memory for the matrix
	matrix<Complex> TijMnTkdagger (dimAfter,dimAfter,0);
	size_t jj;
	size_t kk;
	for ( compressed_matrix<double>::const_iterator1 it1 = Tk.begin1();
	it1 != Tk.end1(); it1++ ) {
		for ( compressed_matrix<double>::const_iterator2 it2 = it1.begin();
		it2 != it1.end(); it2++ ) {
			jj = it2.index2();//transposing Tk jj is index2(), not index1()
			kk = it2.index1();//transposing Tk
			for ( size_t ii = 0; ii < dimAfter; ii++ ) {
	// *it2 is Tk(kk,jj) = trans(Tk)(jj,kk)
	TijMnTkdagger(ii,kk) += TijMn(ii,jj)*(*it2);
			}
		}
	}
	return TijMnTkdagger;
}

vector<Lorentz5Momentum> DensityOperator::boostToRestFrame(const vector<Lorentz5Momentum>& momenta) {
	// We need 2 initial particles
	assert(momenta.size() >= 2);
	// The boosted vectors
	vector<Lorentz5Momentum> vboosted = momenta;

	// The boost should be to the rest frame of the initial particles
	Boost b = (momenta[0] + momenta[1]).findBoostToCM();

	// Boost all of the vectors
	for ( size_t i = 0; i < momenta.size(); i++ ) {
		vboosted[i].boost(b);
	}
		
	return vboosted;
}

bool DensityOperator::compareMomentum(const Lorentz5Momentum& p, const Lorentz5Momentum& q) {
	bool equal = true;
	// Compares two momentum vectors p and q, if they are close enough (defined by the
	// double eps below) it returns true. This is sufficient to distinguish particles
	// in the matrix element as we are not interested in the matrix element for extremely
	// collinear radiation (better described by the parton shower).

	// Create the difference of the two vectors
	const Lorentz5Momentum l = p - q;
	// A relevant size that is guaranteed to be larger than 0
	const Energy2 e2 = p.t()*p.t();
	// Size of the difference that would be considered equal
	const double eps = pow(10.,-15.);

	if ( l.x()*l.x()/e2 > eps )
		equal = false;
	if ( l.y()*l.y()/e2 > eps )
		equal = false;
	if ( l.z()*l.z()/e2 > eps )
		equal = false;
	if ( l.t()*l.t()/e2 > eps )
		equal = false;
	
	return equal;
}




IBPtr DensityOperator::clone() const {
	return new_ptr(*this);
}

IBPtr DensityOperator::fullclone() const {
	return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DensityOperator::persistentOutput(PersistentOStream &) const {
	// *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DensityOperator::persistentInput(PersistentIStream &, int) {
	// *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DensityOperator,HandlerBase>
describeHerwigDensityOperator("Herwig::DensityOperator", "DensityOperator.so");

void DensityOperator::Init() {

	static ClassDocumentation<DensityOperator> documentation
	        ("There is no documentation for the DensityOperator class");

}

