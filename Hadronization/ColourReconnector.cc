// -*- C++ -*-
//
// ColourReconnector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourReconnector class.
//

#include "ColourReconnector.h"
using namespace Herwig;

using CluVecIt = ColourReconnector::CluVecIt;
using Constants::pi;
using Constants::twopi;


DescribeClass<ColourReconnector,Interfaced>
describeColourReconnector("Herwig::ColourReconnector","Herwig.so");

IBPtr ColourReconnector::clone() const {
  return new_ptr(*this);
}

IBPtr ColourReconnector::fullclone() const {
  return new_ptr(*this);
}

void ColourReconnector::rearrange(ClusterVector & clusters) {
  if (_clreco == 0) return;

  // need at least two clusters
  if (clusters.size() < 2) return;
	for (unsigned int i = 0; i < _crIterations; i++) {
		// do the colour reconnection
		switch (_algorithm) {
			case 0:
				_doRecoPlain(clusters);
				break;
			case 1:
				// TODO This algorithm has no dynamic CR
				_doRecoStatistical(clusters);
				break;
			case 2:
				_doRecoBaryonic(clusters);
				break;
			case 3:
				// TODO This algorithm has no dynamic CR
				_doRecoBaryonicMesonic(clusters);
				break;
			case 4:
				_doRecoBaryonicDiquarkCluster(clusters);
				break;
			default:
				assert(false);
		}
	}
}

namespace {
  inline int hasDiquark(const ClusterPtr & cit) {
    int res = 0;
    for (unsigned int i = 0; i<(cit)->numComponents(); i++) {
      if (DiquarkMatcher::Check(*((cit)->particle(i)->dataPtr()))) {
        res++;
      }
    }
    return res;
  }

  // Calculate the rapidity in Lorentz invariant form
  // with respect to n  and nbar (1, -\vec{n}) frame
  // where n = (1, \vec{n}) and nbar (1, -\vec{n}) in COM
  double calculateRelativeRapidity(
      const LorentzVector<double> & n,
      const LorentzVector<double> & nbar,
      const Lorentz5Momentum & q) {
    double rapNum = nbar*q/GeV;
    double rapDen = n*q/GeV;
    double rap = 0.5*(log(rapNum)-log(rapDen));
    return rap;
  }

}


Energy2 ColourReconnector::_clusterMassSum(const PVector & q,
                                           const PVector & aq) const {
  const size_t nclusters = q.size();
  assert (aq.size() == nclusters);
  Energy2 sum = ZERO;
  for (size_t i = 0; i < nclusters; i++)
    sum += ( q[i]->momentum() + aq[i]->momentum() ).m2();
  return sum;
}

double ColourReconnector::_displacement(tcPPtr p, tcPPtr q) const {
  double deltaRap = (p->rapidity() - q->rapidity());
  double deltaPhi = fabs(p->momentum().phi() - q->momentum().phi());
  // keep deltaPhi's below Pi due to periodicity
  if (deltaPhi > M_PI) deltaPhi-=M_PI;
  return sqrt(deltaRap * deltaRap + deltaPhi * deltaPhi);
}

/**
 * Computes circular Mean of three angles alpha, beta, gamma
 * */
static double circularMean(double alpha, double beta, double gamma) {
  double xMean=(cos(alpha)+cos(beta)+cos(gamma))/3.0;
  double yMean=(sin(alpha)+sin(beta)+sin(gamma))/3.0;
  // to make the function fail-save
  if (xMean==0 && yMean==0) return M_PI;
  return atan2(yMean,xMean);
}

namespace {
// ColourFlow for 3 colour flows for baryon state
// ordered permutations by one transposition
// sign of permutations # of difference
// |1> = |123>	+	0
// |2> = |132>	-	1
// |3> = |213>	-	1
// |4> = |231>	+	2
// |5> = |312>	+	2
// |6> = |321>	-	1
int signPermutationState(int i);
int signPermutationState(int i)
{
	if (i==0 || i==3 || i==4) return 1;
	else if (i==1 || i==2 || i==5) return -1;
	else assert(false);
}
// ColourFlow scalar product matrix for 3 
// colour flows in the following basis
// |1> = |123>	+	0
// |2> = |132>	-	1
// |3> = |213>	-	1
// |4> = |231>	+	2
// |5> = |312>	+	2
// |6> = |321>	-	1
unsigned int scalarProducts(int i, int j);
unsigned int scalarProducts(int i, int j)
{
	// Verified for i,j < 3! and 3 colour flows
	if (i>j) return scalarProducts(j,i);
	unsigned int Nc=3;
	if (i==j) return Nc*Nc*Nc;
	switch(i)
	{
		case 0:
		{
			if (j==1 || j==2 || j==5) return Nc*Nc;
			else if (j==3 || j==4) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		case 1:
		{
			if (j==0 || j==3 || j==4) return Nc*Nc;
			else if (j==2 || j==5) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		case 2:
		{
			if (j==0 || j==3 || j==4) return Nc*Nc;
			else if (j==1 || j==5) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		case 3:
		{
			if (j==1 || j==2 || j==5) return Nc*Nc;
			else if (j==0 || j==4) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		case 4:
		{
			if (j==1 || j==2 || j==5) return Nc*Nc;
			else if (j==0 || j==3) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		case 5:
		{
			if (j==0 || j==3 || j==4) return Nc*Nc;
			else if (j==1 || j==2) return Nc;
			else return Nc*Nc*Nc;
			break;
		}
		default:
		assert(false);
	}

	return Nc;
}
}

std::unordered_map<int,double> ColourReconnector::_reconnectionAmplitudesCF2(
    const ClusterPtr & c1, const ClusterPtr & c2, const int) const{
	// Verified according to convention of analytics/matrices2_BCR.nb and not
	// according to https://arxiv.org/abs/1808.06770
	// The same convention as in https://arxiv.org/abs/1808.06770 can be obtained by
	// the mapping 1->1;2->4;3->2;4->3 (here->paper)

	Lorentz5Momentum p1 = c1->colParticle()->momentum();
	Lorentz5Momentum p2 = c1->antiColParticle()->momentum();
	Lorentz5Momentum p3 = c2->colParticle()->momentum();
	Lorentz5Momentum p4 = c2->antiColParticle()->momentum();
	Energy mLightestConstMass = 1000*GeV;
	for (const auto & id : _hadronSpectrum->lightHadronizingQuarks() ) {
		if (getParticleData(id)->constituentMass()<mLightestConstMass)
			mLightestConstMass = getParticleData(id)->constituentMass();
	}
	
	// Scale choice mu = 2*mLightestConstMass*_dynamicCRscale
	double M12 = (p1 + p2).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
	double M24 = (p2 + p4).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
	double M13 = (p1 + p3).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
	double M23 = (p2 + p3).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
	double M14 = (p1 + p4).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
	double M34 = (p3 + p4).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);

	if (
			M12 < 1.0 ||
			M24 < 1.0 ||
			M13 < 1.0 ||
			M23 < 1.0 ||
			M14 < 1.0 ||
			M34 < 1.0
			) {
		throw Exception()
			<< "DynamicCR scale "<< _dynamicCRscale << " too high in ColourReconnector"
			<< " must be less than 1"
			<< " invariant masses:\n"
			<< "M12 = " << sqrt(M12*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< "M24 = " << sqrt(M24*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< "M13 = " << sqrt(M13*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< "M23 = " << sqrt(M23*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< "M14 = " << sqrt(M14*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< "M34 = " << sqrt(M34*sqr(2*mLightestConstMass*_dynamicCRscale))/GeV << "\n"
			<< Exception::eventerror;
	}

	double alphaQCD = _dynamicCRalphaS;

	double logSqrOmega12 = alphaQCD*sqr(log(M12))/(4.0*twopi);
	double logSqrOmega24 = alphaQCD*sqr(log(M24))/(4.0*twopi);
	double logSqrOmega13 = alphaQCD*sqr(log(M13))/(4.0*twopi);
	double logSqrOmega23 = alphaQCD*sqr(log(M23))/(4.0*twopi);
	double logSqrOmega14 = alphaQCD*sqr(log(M14))/(4.0*twopi);
	double logSqrOmega34 = alphaQCD*sqr(log(M34))/(4.0*twopi);

	double Nc = 3.0;
	double U11,U21; // relevant matrix elements
	switch (_dynamicCR)
	{
		case 1:
			{	
				double a = (logSqrOmega34 + logSqrOmega12);
				double b = (logSqrOmega14 + logSqrOmega23);
				double c = (logSqrOmega13 + logSqrOmega24);

				double sqrtDelta = sqrt(Nc*Nc*a*a-4*c*(a+b)-(2*Nc*Nc-4)*a*b+Nc*Nc*b*b+4*c*c);
				U11 = sqrtDelta/tanh(sqrtDelta/2.0)+Nc*(b-a);
				U21 = 2*(c-b);
				break;
			}
		case 2:
			{	
				// not Exponentiated soft anomalous dimension
				// i.e. eq (5.2) Omega matrix in https://arxiv.org/pdf/1808.06770
				U11 = -Nc*(logSqrOmega34+logSqrOmega12);
				U21 =     (logSqrOmega13+logSqrOmega24-(logSqrOmega14+logSqrOmega23));
				break;
			}
		default:
			assert(false);
	}
	
	double TransAmpNoCR      = Nc*(Nc * U11 +      U21);
	double TransAmpMesonicCR = Nc*(U11      + Nc * U21);
	std::unordered_map<int,double> amplitudes;
	// No Colour Reconnection amplitude
	amplitudes[123] = TransAmpNoCR;
	// Mesonic Colour Reconnection amplitude
	amplitudes[213] = TransAmpMesonicCR;
	return amplitudes;
}

/*
 * return the soft gluon evolution probabilities for a two cluster system
 * @returns unitarized probabilities for No, Mesonic and Diquark CR
 * */
std::tuple<double,double,double> ColourReconnector::_dynamicRecoProbabilitiesCF2(
    const ClusterPtr & c1, const ClusterPtr & c2, bool diquarkCR) const{
	static int Nc = 3;
	std::unordered_map<int,double> amplitudes = _reconnectionAmplitudesCF2 (c1, c2); 	
	double pNoCR;
	double pMesonicCR;
	double pDiquarkCR = 0.0;

	double TransAmpNoCR      = sqr(amplitudes[123]);
	double TransAmpMesonicCR = _becomesColour8Cluster(c1,c2) ? 0.0:sqr(amplitudes[213]);

	double sum = TransAmpNoCR + TransAmpMesonicCR;
	double PhaseSpace = 0.0;
	if (diquarkCR && _canMakeDiquarkCluster(c1,c2,PhaseSpace)
      && !hasDiquark(c1) && !hasDiquark(c2)) {
		// Normalization constant of Diquark states to <D|D> = Nc^2 
		static const double ND = sqrt(2*(Nc-1.0)/Nc);
		double TransAmpDiquarkCR = sqr((amplitudes[123]-amplitudes[213])/ND);
    // add phase space suppression for  diquark CR if desired
		sum += _phaseSpaceDiquarkFission ? PhaseSpace*TransAmpDiquarkCR:TransAmpDiquarkCR;
		pDiquarkCR = TransAmpDiquarkCR/sum;
		if (_phaseSpaceDiquarkFission) pDiquarkCR*=PhaseSpace;
		assert( pDiquarkCR<=1.0 && pDiquarkCR>=0.0);
	}
	pNoCR      = TransAmpNoCR/sum;
	pMesonicCR = TransAmpMesonicCR/sum;
	assert( pNoCR<=1.0 && pNoCR>=0.0);
	assert( pMesonicCR<=1.0 && pMesonicCR>=0.0);
	if (_debug)
	{
		ofstream out("WriteOut/kinematicRecoProbability.dat", std::ios::app);
		out << pNoCR << "\t" << pMesonicCR << "\t" << pDiquarkCR << "\t";
		out.close();
	}
	return {pNoCR, pMesonicCR, pDiquarkCR};
}

/*
 * mapping from index of basis vectors to actual permutations
 * */
int ColourReconnector::_stateToPermutation(const int i) const {
	switch (i)
	{
		case 0:
		  // Baryonic CR
			return 0;
		// Mesonic CR's
		case 1:
			return 123;
		case 2:
			return 132;
		case 3:
			return 213;
		case 4:
			return 231;
		case 5:
			return 312;
		case 6:
			return 321;
		default:
			assert(false);
	}
}

/*
 * @returns a Selector object containing the unitarized transition probabilities
 * NOTE: only implemented for 2 and 3 mesonic cluster systems
 * */
Selector<int> ColourReconnector::_selector(const ClusterVector & clusters,
    bool diquarkCR) const {
	switch (clusters.size())
	{
		case 2:
			return getProbabilities2CF(clusters[0], clusters[1], diquarkCR);
			break;
		case 3:
			{
				if (_dynamicCR) {
					return _selectorCF3(clusters, diquarkCR);
				}
				else
				{
					Selector<int> CRoptions;
					double sum = 0;
					const int NpossibilitiesMCR = 6;
					int state;
					for (int i = 0; i < NpossibilitiesMCR; i++)
					{
						state = _stateToPermutation(i+1);
						if (_isColour8Forbidden(state,clusters))
							continue;
						CRoptions.insert(_precoMesonic, state);
						sum += _precoMesonic;
					}

					if (_canMakeBaryonicCluster(clusters[0], clusters[1], clusters[2])) {
						CRoptions.insert(_precoBaryonic, 0);
						sum += _precoBaryonic;
					}
					if (diquarkCR) {
						const int N=3;

						bool first=false;
						// Here i is the index of the quark and j of the antiquark
						// which will be connected to each other (other partons will be 
						// forming the diquark cluster if kinematically viable)
						double PhaseSpace;
						for (int i = 1; i <= N; i++) {
							for (int j = 1; j <= N; j++) {
								state = -(10*i+j);
								if (_isColour8Forbidden(state,clusters))
									continue;
								if (_canMakeDiquarkCluster(
											clusters[i%3]->colParticle(),     clusters[(i+1)%3]->colParticle(),
											clusters[j%3]->antiColParticle(), clusters[(j+1)%3]->antiColParticle(),PhaseSpace)){
                  // Note if _phaseSpaceDiquarkFission is turned off, PhaseSpace == 1.0
									CRoptions.insert(_precoDiquark*PhaseSpace, state);
									if (first){
										sum+=_precoDiquark*PhaseSpace;
										first=false;
									}
								}
							}
						}
					}
					// no CR probability
					CRoptions.insert((1-sum) > 0 ? (1-sum):0, _stateToPermutation(1));
					return CRoptions;
				}
			}
		default:
			assert(false);
	}
	
}

/*
 * @return true if proposal CR `state` will produce a forbidden Octet
 * combination according to _octetOption
 * */
bool ColourReconnector::_isColour8Forbidden(int state, const ClusterVector & clusters) const{
	if (state>0){
		int c1 = ((state/100) % 4) -1;
		int c2 = ((state-(c1+1)*100)/10 % 4) -1;
		int c3 = ((state-(c1+1)*100-(c2+1)*10) % 4) -1;
		assert(c1>=0 && c1<3);
		assert(c2>=0 && c2<3);
		assert(c3>=0 && c3<3);
		assert(c1!=c2 && c1!=c3 && c2 !=c3);
		assert(state==(c1+1)*100+(c2+1)*10+(c3+1));
		if (c1!=0 && _isColour8(clusters[0]->colParticle(),clusters[c1]->antiColParticle()))
			return true;
		if (c2!=1 && _isColour8(clusters[1]->colParticle(),clusters[c2]->antiColParticle()))
			return true;
		if (c3!=2 && _isColour8(clusters[2]->colParticle(),clusters[c3]->antiColParticle()))
			return true;
	}
	else if (state<0){
		int i = -state/10 -1;
		int j = -state - 10*(i+1) - 1;
		assert(-state==((i+1)*10+(j+1)));
		if (_isColour8(clusters[i]->colParticle(),clusters[j]->antiColParticle()))
			return true;
		
	}
	return false;
}

/*
 * @returns Selector object with unitarized probabilities for 3 mesonic cluster system
 * */
Selector<int> ColourReconnector::_selectorCF3(const ClusterVector & clusters, bool diquarkCR) const {
	std::unordered_map<int, double> amplitudes = _reconnectionAmplitudesSGE(clusters);
	double sum2 = 0;
	double sumBaryon = 0;
	static const int Nc = 3;
	// Normalization for baryonic state such that <B|B> = Nc^3
	static const double NB = sqrt(6*(Nc*Nc-3*Nc+2)/(Nc*Nc));
	Selector<int> CRoptions;
	double amp2;
	const int NpossibilitiesMCR = 6;
	int state;
	bool canMakeBaryon = _canMakeBaryonicCluster(clusters[0], clusters[1], clusters[2]);
	for (int i = 0; i < NpossibilitiesMCR; i++)
	{
		state=_stateToPermutation(i+1);
		if (canMakeBaryon)
			sumBaryon +=     amplitudes[state]*signPermutationState(i)/NB;
		if (_isColour8Forbidden(state,clusters))
			continue;
		amp2 =       sqr(amplitudes[state]);
		sum2 += amp2;
		CRoptions.insert(amp2, state);
	}
	sum2+=sumBaryon*sumBaryon;
	if (diquarkCR) {
		// Normalization constant of Diquark states to <D|D> = Nc^2 
		static const double ND = sqrt(2*(Nc-1.0)/Nc);
		const int N=3;
		double PhaseSpace;

		int perm1;
		int perm2;
		// Here i is the index of the quark and j of the antiquark
		// which will be connected to each other (other partons will be 
		// forming the diquark cluster if kinematically viable)
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				assert(i%3!=(i-1));
				assert((i+1)%3!=(i-1));
				assert(j%3!=(j-1));
				assert((j+1)%3!=(j-1));
				if (_canMakeDiquarkCluster(
						clusters[i%3]->colParticle(),     clusters[(i+1)%3]->colParticle(),
						clusters[j%3]->antiColParticle(), clusters[(j+1)%3]->antiColParticle()	,PhaseSpace))
				{
					state = -(10*i+j);
					if (_isColour8Forbidden(state,clusters))
						continue;
					switch (i)
					{
						case 1:
							{
								perm1 = 100*j + 10*(((j+1)%3)+1) + (((j)%3)+1);
								perm2 = 100*j + 10*(((j)%3)+1)   + (((j+1)%3)+1);
								amp2 = sqr((amplitudes[perm1] - amplitudes[perm2])/ND);
								sum2 += amp2;
                if (amp2==0 || amplitudes[perm2]==0 || amplitudes[perm1]==0 ){
                  throw Exception()
                    << "Zero transition amplitudes in ColourReconnector::_selectorCF3\n"
                    << "amp23 = "<< amp2 << "\n"
                    << "perm13 = "<< perm1 << "\n"
                    << "perm23 = "<< perm2 << "\n"
                    << Exception::warning;
                }
								assert((((j)%3)+1)!=j);
								assert((((j+1)%3)+1)!=j);
								CRoptions.insert(PhaseSpace*amp2, state);
								break;
							}
						case 2:
							{
								perm1 = 100*(((j)%3)+1)   + 10*j + (((j+1)%3)+1);
								perm2 = 100*(((j+1)%3)+1) + 10*j + (((j)%3)+1);
								amp2 = sqr((amplitudes[perm1] - amplitudes[perm2])/ND);
								sum2 += amp2;
								if (amp2==0 || amplitudes[perm2]==0 || amplitudes[perm1]==0 ){
                  throw Exception()
                    << "Zero transition amplitudes in ColourReconnector::_selectorCF3\n"
                    << "amp23 = "<< amp2 << "\n"
                    << "perm13 = "<< perm1 << "\n"
                    << "perm23 = "<< perm2 << "\n"
                    << Exception::warning;
								}
								assert((((j)%3)+1)!=j);
								assert((((j+1)%3)+1)!=j);
								CRoptions.insert(PhaseSpace*amp2, state);
								break;
							}
						case 3:
							{
								perm1 = 100*(((j+1)%3)+1) + 10*(((j)%3)+1)   + j;
								perm2 = 100*(((j)%3)+1)   + 10*(((j+1)%3)+1) + j;
								amp2 = sqr((amplitudes[perm1] - amplitudes[perm2])/ND);
								sum2 += amp2;
								if (amp2==0 || amplitudes[perm2]==0 || amplitudes[perm1]==0 ){
                  throw Exception()
                    << "Zero transition amplitudes in ColourReconnector::_selectorCF3\n"
                    << "amp23 = "<< amp2 << "\n"
                    << "perm13 = "<< perm1 << "\n"
                    << "perm23 = "<< perm2 << "\n"
                    << Exception::warning;
								}
								assert((((j)%3)+1)!=j);
								assert((((j+1)%3)+1)!=j);
								CRoptions.insert(PhaseSpace*amp2, state);
								break;
							}
						default:
							assert(false);
					}
				}
			}
		}
	}
	if (canMakeBaryon)
		CRoptions.insert(sumBaryon*sumBaryon, _stateToPermutation(0));
	return CRoptions;	
}

/*
 * @returns a hash table of the transition amplitudes (NOT the probabilities)
 * */
std::unordered_map<int, double> ColourReconnector::_reconnectionAmplitudesSGE(
    const ClusterVector & clusters) const {
	int size = clusters.size();
	assert(clusters.size()<4);
	switch(size){
		case 2:
			{
				const std::unordered_map<int, double> & amplitudes = _reconnectionAmplitudesCF2(clusters[0],clusters[1]);
				return amplitudes;
			}
		case 3:
			{
				std::unordered_map<int, double> amplitudes;
				if (clusters[0]->numComponents()!=2 ||
						clusters[1]->numComponents()!=2	||
						clusters[2]->numComponents()!=2	||
						hasDiquark(clusters[0]) ||
						hasDiquark(clusters[1]) ||
						hasDiquark(clusters[2])   ) {
					std::cout << "reject config. should reject before somehow\n";
					return amplitudes;
				}
				const int N = 6; // 2*Nclu;
				double omIJ[N][N];
				double alphaQCD=_dynamicCRalphaS;
				Lorentz5Momentum mom_i, mom_j;
				// Ordering of omIJ and scalarProds is according to {clu1_col, clu1_anti, clu2_col, clu2_anti,...}
				Energy mLightestConstMass = 1000*GeV;
				for (auto id : _hadronSpectrum->lightHadronizingQuarks() ) {
					if (getParticleData(id)->constituentMass()<mLightestConstMass)
						mLightestConstMass = getParticleData(id)->constituentMass();
				}
				for (int i = 0; i < N; i++)
				{
					if ((i+1)%2==1) mom_i=clusters[i/2]->colParticle()->momentum();
					else mom_i=clusters[(i-1)/2]->antiColParticle()->momentum();
					for (int j = i+1; j < N; j++)
					{
						if ((j+1)%2==1) mom_j=clusters[j/2]->colParticle()->momentum();
						else mom_j=clusters[(j-1)/2]->antiColParticle()->momentum();
						double rho = (mom_i+mom_j).m2()/sqr(2*mLightestConstMass*_dynamicCRscale);
						if (rho < 1.0) {
							throw Exception()
								<< "DynamicCR scale "<< _dynamicCRscale << " too high in ColourReconnector"
								<< " must be less than 1"
								<< " Found parton invariant mass combinations with invariant mass ratio "<< (mom_i+mom_j).m2()/sqr((mom_i.m()+mom_j.m()))
								<< Exception::eventerror;
						}
						omIJ[i][j]=alphaQCD*sqr(log(rho))/(4.0*twopi);
					}
				}
				boost::numeric::ublas::matrix<double> * Uevolve = new boost::numeric::ublas::matrix<double>(N,N);
				boost::numeric::ublas::matrix<double> Omega(N,N);
				int Nc = 3;
				// Verified Omega Matrix with analytics/matrices2_BCR.nb 
				Omega(0,0) = - Nc * (omIJ[0][1]+omIJ[2][3]+omIJ[4][5]);
				Omega(1,0) =   (omIJ[2][4]-omIJ[3][4]-omIJ[2][5]+omIJ[3][5]);
				Omega(2,0) =   (omIJ[0][2]-omIJ[1][2]-omIJ[0][3]+omIJ[1][3]);
				Omega(3,0) = 0.0;
				Omega(4,0) = 0.0;
				Omega(5,0) =   (omIJ[0][4]-omIJ[1][4]-omIJ[0][5]+omIJ[1][5]);

				Omega(0,1) = - (omIJ[2][3]-omIJ[2][4]-omIJ[3][5]+omIJ[4][5]);
				Omega(1,1) = - Nc * (omIJ[0][1]+omIJ[3][4]+omIJ[2][5]);
				Omega(2,1) = 0.0;
				Omega(3,1) = - (omIJ[0][3]-omIJ[1][3]-omIJ[0][4]+omIJ[1][4]);
				Omega(4,1) =   (omIJ[0][2]-omIJ[1][2]-omIJ[0][5]+omIJ[1][5]);
				Omega(5,1) = 0.0;

				Omega(0,2) = - (omIJ[0][1]-omIJ[0][2]-omIJ[1][3]+omIJ[2][3]);
				Omega(1,2) = 0.0;
				Omega(2,2) = - Nc * (omIJ[1][2]+omIJ[0][3]+omIJ[4][5]);
				Omega(3,2) = - (omIJ[1][4]-omIJ[2][4]-omIJ[1][5]+omIJ[2][5]);
				Omega(4,2) =   (omIJ[0][4]-omIJ[3][4]-omIJ[0][5]+omIJ[3][5]);
				Omega(5,2) = 0.0;

				Omega(0,3) = 0.0;
				Omega(1,3) = - (omIJ[0][1]-omIJ[1][3]-omIJ[0][4]+omIJ[3][4]);
				Omega(2,3) = - (omIJ[1][2]-omIJ[2][4]-omIJ[1][5]+omIJ[4][5]);
				Omega(3,3) = - Nc * (omIJ[0][3]+omIJ[1][4]+omIJ[2][5]);
				Omega(4,3) = 0.0;
				Omega(5,3) =   (omIJ[0][2]-omIJ[2][3]-omIJ[0][5]+omIJ[3][5]);

				Omega(0,4) = 0.0;
				Omega(1,4) = - (omIJ[0][1]-omIJ[0][2]-omIJ[1][5]+omIJ[2][5]);
				Omega(2,4) = - (omIJ[0][3]-omIJ[0][4]-omIJ[3][5]+omIJ[4][5]);
				Omega(3,4) = 0.0;
				Omega(4,4) = - Nc * (omIJ[1][2]+omIJ[3][4]+omIJ[0][5]);
				Omega(5,4) =   (omIJ[1][3]-omIJ[2][3]-omIJ[1][4]+omIJ[2][4]);

				Omega(0,5) = - (omIJ[0][1]-omIJ[0][4]-omIJ[1][5]+omIJ[4][5]);
				Omega(1,5) = 0.0;
				Omega(2,5) = 0.0;
				Omega(3,5) =   (omIJ[0][2]-omIJ[0][3]-omIJ[2][5]+omIJ[3][5]);
				Omega(4,5) = - (omIJ[1][2]-omIJ[1][3]-omIJ[2][4]+omIJ[3][4]);
				Omega(5,5) = - Nc * (omIJ[2][3]+omIJ[1][4]+omIJ[0][5]);
				switch (_dynamicCR)
				{
					case 1:
						// Exponentiated
						*Uevolve = expm_pad(Omega);
						break;
					case 2:
						// NotExponentiated
						*Uevolve = Omega;
						break;
					default:
						assert(false);
				}
				// std::cout << *Uevolve << std::endl;
				double amp1toJ;
				for (int J = 0; J < N; J++)
				{
					// amp1toJ is here for each J transition amplitude 1->J
					amp1toJ=0;

					for (int i = 0; i < N; i++)
					{
						// amp1toJ corresponds to the Uevolve operator applied to the |1> state
						// and projected by fixed <J|
						// The resulting amp1toJ is <J|U|1>
						// TODO : Generalize here to general starting states i.e. |1> -> |initial>
						amp1toJ+=(*Uevolve)(i,0)*scalarProducts(i,J);
					}
					if (std::isnan(amp1toJ) || std::isinf(amp1toJ)){
						throw Exception() << "nan or inf transition probability in ColourReconnector::_reconnectionAmplitudesSGE"
							<< Exception::runerror;
					}
					amplitudes[_stateToPermutation(J+1)] = amp1toJ;
				}
				delete Uevolve;
				return amplitudes;
			}
		default:
			{
				throw Exception() << "Found cluster set of " << size << " in ColourReconnector::_reconnectionAmplitudesSGE (can only handle 2 or 3 colour Flows)"
				<< Exception::runerror;
			}
	}
	return std::unordered_map<int, double>();
}

double ColourReconnector::_displacementBaryonic(tcPPtr q1, tcPPtr q2, tcPPtr q3) const {
  if (_junctionMBCR) {
    /**
     * Junction-like option i.e. displacement
     * from "junction centre" (mean rapidity/phi)
     */
    double rap1=q1->rapidity();
    double rap2=q2->rapidity();
    double rap3=q3->rapidity();

    double phi1=q1->momentum().phi();
    double phi2=q2->momentum().phi();
    double phi3=q3->momentum().phi();
    double meanRap=(rap1 + rap2 + rap3)/3.0;
	// Use circularMean for defining a sensible mean of an angle
	double meanPhi=circularMean(phi1,phi2,phi3);

	double deltaPhi1=fabs(phi1-meanPhi);
	double deltaPhi2=fabs(phi2-meanPhi);
	double deltaPhi3=fabs(phi3-meanPhi);
	// keep deltaPhi's below Pi due to periodicity
	if (deltaPhi1>M_PI) deltaPhi1-=M_PI;
	if (deltaPhi2>M_PI) deltaPhi2-=M_PI;
	if (deltaPhi3>M_PI) deltaPhi3-=M_PI;
	double delR;

	delR  = sqrt( (rap1-meanRap)*(rap1-meanRap) + deltaPhi1*deltaPhi1 );
	delR += sqrt( (rap2-meanRap)*(rap2-meanRap) + deltaPhi2*deltaPhi2 );
	delR += sqrt( (rap3-meanRap)*(rap3-meanRap) + deltaPhi3*deltaPhi3 );
    return delR;
  } else {
    /* just summing up all possible 2 quark displacements */
    return _displacement(q1, q2) + _displacement(q1, q3) + _displacement(q2, q3);
  }
}
bool ColourReconnector::_containsColour8(const ClusterVector & cv,
                                         const vector<size_t> & P) const {
  assert (P.size() == cv.size());
  for (size_t i = 0; i < cv.size(); i++) {
    tcPPtr p = cv[i]->colParticle();
    tcPPtr q = cv[P[i]]->antiColParticle();
    if (_isColour8(p, q)) return true;
  }
  return false;
}


void ColourReconnector::_doRecoStatistical(ClusterVector & cv) const {

  const size_t nclusters = cv.size();

  // initially, enumerate (anti)quarks as given in the cluster vector
  ParticleVector q, aq;
	q.reserve(nclusters);
	aq.reserve(nclusters);
  for (size_t i = 0; i < nclusters; i++) {
    q.push_back( cv[i]->colParticle() );
    aq.push_back( cv[i]->antiColParticle() );
  }

  // annealing scheme
  Energy2 t, delta;
  Energy2 lambda = _clusterMassSum(q,aq);
  const unsigned _ntries = _triesPerStepFactor * nclusters;

  // find appropriate starting temperature by measuring the largest lambda
  // difference in some dry-run random rearrangements
  {
    vector<Energy2> typical;
    typical.reserve(10);
    for (int i = 0; i < 10; i++) {
      const pair <int,int> toswap = _shuffle(q,aq,5);
      ParticleVector newaq = aq;
      swap (newaq[toswap.first], newaq[toswap.second]);
      Energy2 newlambda = _clusterMassSum(q,newaq);
      typical.push_back( abs(newlambda - lambda) );
    }
    t = _initTemp * Math::median(typical);
  }

  // anneal in up to _annealingSteps temperature steps
  for (unsigned step = 0; step < _annealingSteps; step++) {

    // For this temperature step, try to reconnect _ntries times. Stop the
    // algorithm if no successful reconnection happens.
    unsigned nSuccess = 0;
    for (unsigned it = 0; it < _ntries; it++) {

      // make a random rearrangement
      const unsigned maxtries = 10;
      const pair <int,int> toswap = _shuffle(q,aq,maxtries);
      const int i = toswap.first;
      const int j = toswap.second;

      // stop here if we cannot find any allowed reconfiguration
      if (i == -1) break;

      // create a new antiquark vector with the two partons swapped
      ParticleVector newaq = aq;
      swap (newaq[i], newaq[j]);

      // Check if lambda would decrease. If yes, accept the reconnection. If no,
      // accept it only with a probability given by the current Boltzmann
      // factor. In the latter case we set p = 0 if the temperature is close to
      // 0, to avoid division by 0.
      Energy2 newlambda = _clusterMassSum(q,newaq);
      delta = newlambda - lambda;
      double prob = 1.0;
      if (delta > ZERO) prob = ( abs(t) < 1e-8*MeV2 ) ? 0.0 : exp(-delta/t);
      if (UseRandom::rnd() < prob) {
        lambda = newlambda;
        swap (newaq, aq);
        nSuccess++;
      }
    }
    if (nSuccess == 0) break;

    // reduce temperature
    t *= _annealingFactor;
  }

  // construct the new cluster vector
  ClusterVector newclusters;
	newclusters.reserve(nclusters);
  for (size_t i = 0; i < nclusters; i++) {
    ClusterPtr cl = new_ptr( Cluster( q[i], aq[i] ) );
    newclusters.push_back(cl);
  }
  swap(newclusters,cv);
  return;
}

/*
 * @returns Selector with unitarized CR probabilities for a 2 mesonic cluster system
 * */
Selector <int> ColourReconnector::getProbabilities2CF(const ClusterPtr & c1, const ClusterPtr & c2, bool diquarkCR) const{
	Selector <int> res;
	if (_dynamicCR) {
		double pNoCR;
		double pMCR;
		double pDCR;
		std::tie(pNoCR, pMCR, pDCR) = _dynamicRecoProbabilitiesCF2(c1, c2, diquarkCR);
		// Mesonic CR
		if ( !(   _isColour8(c1->colParticle(), c2->antiColParticle())
				   || _isColour8(c2->colParticle(), c1->antiColParticle()) ) ) {
			res.insert(pMCR, 213);
		}
		if (diquarkCR && pDCR>0.0  && _canMakeDiquarkCluster(c1,c2)) {
			// Diquark CR
			res.insert(pDCR, -213);
		}
		// No CR
		res.insert(pNoCR, 123);
	}
	else {
		double wMCR = _precoMesonic;
		double sum = 0.0;
		// Mesonic CR
		if ( !(   _isColour8(c1->colParticle(), c2->antiColParticle())
				   || _isColour8(c2->colParticle(), c1->antiColParticle()) ) ) {
			res.insert(wMCR, 213);
			sum+=wMCR;
		}
		double PhaseSpace;
		if (diquarkCR && _canMakeDiquarkCluster(c1,c2,PhaseSpace)) {
			double wDCR = _precoDiquark*PhaseSpace;
			// Diquark CR
			res.insert(wDCR, -213);
			sum+=wDCR;
		}
		// No CR
		res.insert((1.0-sum)>0 ? (1.0-sum):0.0, 123);
	}
	return res;
}

void ColourReconnector::_doRecoPlain(ClusterVector & cv) const {

  ClusterVector newcv = cv;

	ClusterVector deleted; deleted.reserve(cv.size());
  // try to avoid systematic errors by randomising the reconnection order
  long (*p_irnd)(long) = UseRandom::irnd;
  random_shuffle( newcv.begin(), newcv.end(), p_irnd );

  // iterate over all clusters
  for (CluVecIt cit = newcv.begin(); cit != newcv.end(); cit++) {
    // find the cluster which, if reconnected with *cit, would result in the
    // smallest sum of cluster masses
		// skip the diquark clusters to be deleted later 2->1 cluster
		if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
			continue;
		// skip diquark clusters
    if ((*cit)->numComponents()>2 || (_diquarkCR && hasDiquark(*cit))) continue;
    // NB this method returns *cit if no reconnection partner can be found
    CluVecIt candidate = _dynamicCR ? _findRecoPartnerPlainDynamic(cit, newcv, deleted, _diquarkCR>0):_findRecoPartnerPlain(cit, newcv, deleted);


    // skip this cluster if no possible reshuffling partner can be found
    if (candidate == cit) continue;

    // accept the reconnection with probability PrecoProb.
		// ClusterVector cluvec = {*cit, *candidate};
		// const Selector <int> & sel = getProbabilities2CF(*cit, *candidate, _diquarkCR>0);
		// const Selector <int> & sel = _selector(cluvec, _diquarkCR>0);
		const Selector <int> & sel = _selector({*cit, *candidate}, _diquarkCR>0);
		enum Selection2ColourFlows {
			NoReconnection = 123,
			MesonicReconnection = 213,
			DiquarkReconnection = -213,
		};
		switch (sel.select(UseRandom::rnd()))
		{
			case MesonicReconnection:
				{
					// Mesonic Colour Reconnection
					pair <ClusterPtr,ClusterPtr> reconnected = _reconnect(*cit, *candidate);

          *cit = reconnected.first;
          *candidate = reconnected.second;
					break;
				}
			case DiquarkReconnection:
				{
					ClusterPtr DiqCluster;
					if (_makeDiquarkCluster(*cit, *candidate, DiqCluster)){
						deleted.push_back(*candidate);

						// Note that these must be the cit,candidate and not cluvec
						*cit = DiqCluster;
					}
					break;
				}
			case NoReconnection:
				// No colour Reconnection
				break;
			default:
				assert(false);
		}

  }

	if (deleted.size()==0) {
		cv.swap(newcv);
	}
	else {
		// create a new vector of clusters except for the ones which are "deleted" during
		// Diquark reconnection
		ClusterVector clustervector;
		clustervector.reserve(newcv.size());
		for (const auto & cluster : newcv)
			if (find(deleted.begin(),
						deleted.end(), cluster) == deleted.end())
				clustervector.push_back(cluster);

		cv.swap(clustervector);
	}
}



void ColourReconnector::_doRecoBaryonicDiquarkCluster(ClusterVector & cv) const {
	ClusterVector newcv = cv;

	ClusterVector deleted; deleted.reserve(cv.size());

	// try to avoid systematic errors by randomising the reconnection order
	long (*p_irnd)(long) = UseRandom::irnd;
	random_shuffle(newcv.begin(), newcv.end(), p_irnd);

	Selector <int> sel;
	ClusterVector cluvec;
	cluvec.reserve(3);
	int selection;
	// iterate over all clusters
	for (CluVecIt cit = newcv.begin(); cit != newcv.end(); ++cit) {
		//avoid clusters already containing diuarks
		if (hasDiquark(*cit)) continue;

		//skip the cluster to be deleted later 3->2 cluster
		if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
			continue;

		// Skip all found baryonic and Tetra clusters, this biases the 
		// algorithm but implementing something like re-reconnection
		// is ongoing work
		if ((*cit)->numComponents()>=3) continue;

		// Find a candidate suitable for reconnection
		CluVecIt candidate1, candidate2;
		unsigned typeOfReconnection = 0;
    // rapidity based similar to Baryonic
    _findPartnerBaryonicDiquarkCluster(cit, newcv,
        typeOfReconnection,
        deleted,
        candidate1,
        candidate2);

		switch (typeOfReconnection)
		{
			case 0:
				// No CR found
				continue;
			case 3:
			case 1:
				// Mesonic or Diquark CR with 2 Colour Flows
				cluvec = {*cit, *candidate1};
				break;
			case 2:
				// Baryonic CR with 3 Colour Flows
				cluvec = {*cit, *candidate1, *candidate2};
				break;
			default:
				assert(false);
		}
		if (candidate2 != cit) {
			cluvec = {*cit, *candidate1, *candidate2};
		}
		else if (candidate1 != cit) {
			cluvec = {*cit, *candidate1};
		}
		else {
		  // no CR
			continue;
		}
		if (_dynamicCR) {
			sel = _selector(cluvec, _diquarkCR>0);
			if (typeOfReconnection!=0 && sel.empty()){
				throw Exception()
					<< "No Selection availible"
					<< "ColourReconnector::_doRecoBaryonicDiquarkCluster()"
					<< Exception::runerror;
			}
			selection = sel.select(UseRandom::rnd());
			sel.clear();
			assert(sel.empty());
		}
		else {
			switch (typeOfReconnection)
			{
				case 0:
					// No CR
					std::cout << "Should never execute!!!" << std::endl;
					selection = 123;
					break;
				case 1:
					// Mesonic CR with 2 Colour Flows
					if (UseRandom::rnd() < _precoMesonic)
						selection = 213;
					else
						selection = 123;
					break;
				case 2:
					// Baryonic CR with 3 Colour Flows
					if (UseRandom::rnd() < _precoBaryonic
							&& _canMakeBaryonicCluster(*cit,*candidate1,*candidate2))
						selection = 0;
					else
						selection = 123;
					break;
				case 3:
					// Diquark CR with 2 Colour Flows
					if (UseRandom::rnd() < _precoDiquark)
						selection = -1000;
					else
						selection = 123;
					break;
				default:
					assert(false);
			}
		}
		if (selection == 0){
			// Baryonic CR (only one option for 3 colourflows)
			deleted.push_back(*candidate2);
			// Function that does the reconnection from 3 -> 2 clusters
			ClusterPtr b1, b2;
			_makeBaryonicClusters(*cit,*candidate1,*candidate2, b1, b2);
			// Note that these must be the cit,candidate1 and not cluvec
			*cit = b1;
			*candidate1 = b2;
		}
		else if (selection > 0) {
			// Mesonic CR
			int c1 = ((selection/100) % 4) -1;
			int c2 = ((selection-(c1+1)*100)/10 % 4) -1;
			int c3 = ((selection-(c1+1)*100-(c2+1)*10) % 4) -1;
			if (   c1==0
					&& c2==1
					&& c3==2){
				// noCR
				continue;
			}
			assert(c1>=0 && c1<3);
			assert(c2>=0 && c2<3);
			assert(c3>=0 && c3<3);
			assert(c1!=c2 && c1!=c3 && c2 !=c3);
			// if last cluster is untouched we do only reconnect the first two clusters
			if (c3==2) {
				const auto & reconnected = _reconnect(*cit,*candidate1);
				// Note that these must be the cit, candidate1 and not cluvec
				*cit = reconnected.first;
				*candidate1 = reconnected.second;
			}
			else {
				// Form clusters (0,c1) (1,c2) (2,c3)
				const int infoMCR[3] = {c1, c2, c3}; 
				const auto & reconnected = _reconnect3Mto3M(*cit,*candidate1,*candidate2, infoMCR);

				// Note that these must be the cit,candidateX and not cluvec
				*cit = std::get<0>(reconnected);
				*candidate1 = std::get<1>(reconnected);
				*candidate2 = std::get<2>(reconnected);
			}
		}
		else if (selection < 0) {
			// We will delete the candidate1 mesonic clusters
			// need to check with 2CF solution
			// to form a diquark cluster
			if (cluvec.size()==3 && -selection<999) {
				const auto & reconnected = _reconnect3MtoMD(cluvec, -selection);
				*cit = std::get<0>(reconnected);
				*candidate1 = std::get<1>(reconnected);
				deleted.push_back(*candidate2);
			}
			else if (cluvec.size()==2 || -selection>999){
				ClusterPtr DiqCluster;
				if (_makeDiquarkCluster(*cit, *candidate1, DiqCluster)){
					deleted.push_back(*candidate1);

					// Note that these must be the cit,candidate and not cluvec
					*cit = DiqCluster;
				}
			}
			else {
				assert(false);
			}
		}
		else {
			std::cout << "\nError in selection = "<<selection<<"\n" << std::endl;
			assert(false);
		}
	}
	ClusterVector addedcv;

	// create a new vector of clusters except for the ones which are "deleted" during
	// baryonic reconnection
	ClusterVector clustervector;
	clustervector.reserve(addedcv.size()+newcv.size());
	// add new clusters
	for (const auto & c : addedcv)
		clustervector.push_back(c);
	// delete deleted clusters
	for (const auto & cluster : newcv)
		if (find(deleted.begin(),
					deleted.end(), cluster) == deleted.end())
			clustervector.push_back(cluster);

	swap(cv, clustervector);



}

// Implementation of the baryonic reconnection algorithm
void ColourReconnector::_doRecoBaryonic(ClusterVector & cv) const {

  ClusterVector newcv = cv;

  ClusterVector deleted; deleted.reserve(cv.size());

  // try to avoid systematic errors by randomising the reconnection order
  long (*p_irnd)(long) = UseRandom::irnd;
  random_shuffle(newcv.begin(), newcv.end(), p_irnd);

  double ProbabilityMesonic  = _precoMesonic;
  double ProbabilityBaryonic = _precoBaryonic;
	ClusterVector cluvec;
	cluvec.reserve(3);
  // iterate over all clusters
  for (CluVecIt cit = newcv.begin(); cit != newcv.end(); ++cit) {
		//avoid clusters already containing diuarks
		if (hasDiquark(*cit)) continue;

		//skip the cluster to be deleted later 3->2 cluster
		if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
			continue;

		// Skip all found baryonic clusters, this biases the algorithm but implementing
		// something like re-reconnection is ongoing work
		if ((*cit)->numComponents()>=3) continue;

		// Find a candidate suitable for reconnection
		CluVecIt baryonic1, baryonic2;
		bool isBaryonicCandidate = false;
		CluVecIt candidate = _findPartnerBaryonic(cit, newcv,
				isBaryonicCandidate,
				deleted,
				baryonic1, baryonic2);

		// skip this cluster if no possible reconnection partner can be found
		if ( !isBaryonicCandidate && *candidate==*cit)
			continue;
		if (_dynamicCR) {
			cluvec.clear();
			if (isBaryonicCandidate) {
				cluvec={*cit,*baryonic1,*baryonic2};
			}
			else{
				cluvec={*cit,*candidate};
			}
			const Selector <int> & sel = _selector(cluvec, _diquarkCR>0);
			int selection = sel.select(UseRandom::rnd());
			// 3 aligned meson case
			// Normal 2->2 Colour reconnection
			if (selection == 0){
				// Baryonic CR only one option
				deleted.push_back(*baryonic2);
				// Function that does the reconnection from 3 -> 2 clusters
				ClusterPtr b1, b2;
				_makeBaryonicClusters(*cit,*baryonic1,*baryonic2, b1, b2);

				*cit = b1;
				*baryonic1 = b2;
			}
			else if (selection > 0) {
				if (isBaryonicCandidate) {
					int c1 = ((selection/100) % 4) -1;
					int c2 = ((selection-(c1+1)*100)/10 % 4) -1;
					int c3 = ((selection-(c1+1)*100-(c2+1)*10) % 4) -1;
					if (   c1==0
							&& c2==1
							&& c3==2){
						// noCR
						continue;
					}
					assert(c1>=0 && c1<3);
					assert(c2>=0 && c2<3);
					assert(c3>=0 && c3<3);
					assert(c1!=c2 && c1!=c3 && c2 !=c3);
					if (c3==2 ) {
						const auto & reconnected = _reconnect(*cit,*baryonic1);
						*cit = reconnected.first;
						*baryonic1 = reconnected.second;
					}
					else {
						// Form clusters (0,c1) (1,c2) (2,c3)
						const int infoMCR[3] = {c1, c2, c3}; 
						const auto & reconnected = _reconnect3Mto3M(*cit,*baryonic1,*baryonic2, infoMCR);

						*cit = std::get<0>(reconnected);
						*baryonic1 = std::get<1>(reconnected);
						*baryonic2 = std::get<2>(reconnected);
					}
				}
				else {
					const auto & reconnected = _reconnect(*cit,*candidate);
					*cit = reconnected.first;
					*candidate = reconnected.second;
				}
			}
			else if (selection < 0) {
				// We will delete the candidate1 mesonic clusters
				// need to check with 2CF solution
				// to form a diquark cluster

				if (cluvec.size()==3) {
					const auto & reconnected = _reconnect3MtoMD(cluvec, -selection);
					*cit = std::get<0>(reconnected);
					*baryonic1 = std::get<1>(reconnected);
					deleted.push_back(*baryonic2);
				}
				else if (cluvec.size()==2) {
					ClusterPtr DiqCluster;
					if (_makeDiquarkCluster(*cit, *candidate, DiqCluster)){
						deleted.push_back(*candidate);

						// Note that these must be the cit,candidate and not cluvec
						*cit = DiqCluster;
					}
				}
			}
			else {
				std::cout << "\nError in selection = "<<selection<<"\n" << std::endl;
				assert(false);
			}
		}
		else {
			// 3 aligned meson case
			if ( isBaryonicCandidate
					&& UseRandom::rnd() < ProbabilityBaryonic
					&& _canMakeBaryonicCluster(*cit, *baryonic1, *baryonic2)) {
				deleted.push_back(*baryonic2);

				// Function that does the reconnection from 3 -> 2 clusters
				ClusterPtr b1, b2;
				_makeBaryonicClusters(*cit, *baryonic1, *baryonic2, b1, b2);

				*cit = b1;
				*baryonic1 = b2;

				// Baryonic2 is easily skipped in the next loop
			}

			// Normal 2->2 Colour reconnection
			if (!isBaryonicCandidate
					&& UseRandom::rnd() < ProbabilityMesonic) {
				const auto & reconnected = _reconnect(*cit, *candidate);
				*cit = reconnected.first;
				*candidate = reconnected.second;
			}
		}
  }

  // create a new vector of clusters except for the ones which are "deleted" during
  // baryonic reconnection
  ClusterVector clustervector;
  for (const auto & cluster : newcv)
    if (find(deleted.begin(),
             deleted.end(), cluster) == deleted.end())
      clustervector.push_back(cluster);

  swap(cv, clustervector);
}


bool ColourReconnector::_clustersFarApart( const std::vector<CluVecIt> & clu ) const {
  int Ncl=clu.size();
  assert(Ncl<=3);
  if (Ncl==1) {
    return false;
  } else if (Ncl==2) {
	// veto if Clusters further apart than _maxDistance
	if (_localCR && ((*clu[0])->vertex().vect()-(*clu[1])->vertex().vect()).mag() > _maxDistance) return true;
	// veto if Clusters have negative spacetime difference
	if (_causalCR && ((*clu[0])->vertex()-(*clu[1])->vertex()).m() < ZERO) return true;
  } else if (Ncl==3) {
	// veto if Clusters further apart than _maxDistance
	if (_localCR && ((*clu[0])->vertex().vect()-(*clu[1])->vertex().vect()).mag() > _maxDistance) return true;
	if (_localCR && ((*clu[1])->vertex().vect()-(*clu[2])->vertex().vect()).mag() > _maxDistance) return true;
	if (_localCR && ((*clu[0])->vertex().vect()-(*clu[2])->vertex().vect()).mag() > _maxDistance) return true;
	// veto if Clusters have negative spacetime difference
	if (_causalCR && ((*clu[0])->vertex()-(*clu[1])->vertex()).m() < ZERO) return true;
	if (_causalCR && ((*clu[1])->vertex()-(*clu[2])->vertex()).m() < ZERO) return true;
	if (_causalCR && ((*clu[0])->vertex()-(*clu[2])->vertex()).m() < ZERO) return true;
  }

  return false;
}



void ColourReconnector::_doReco2BeamClusters(ClusterVector & cv) const {
  // try other option
  tPPtr p1Di=(cv[0])->colParticle();
  tPPtr p2Di=(cv[1])->colParticle();

  tPPtr p1Q=(cv[0])->antiColParticle();
  tPPtr p2Q=(cv[1])->antiColParticle();

  double min_dist=_displacement(p1Di,p1Q)+_displacement(p2Di,p2Q);

  if ((_displacement(p1Di,p2Q)+_displacement(p1Di,p2Q))<min_dist) {
    _reconnect(cv[0],cv[1]);
  }
  return;
}



void ColourReconnector::_doRecoBaryonicMesonic(ClusterVector & cv) const {
  if (cv.size() < 3) {
    return;
  }

  ClusterVector newcv = cv;
  newcv.reserve(2*cv.size());

  ClusterVector deleted;
  deleted.reserve(cv.size());

  // counters for numbers of mesons and baryons selected
  unsigned num_meson = 0;
  unsigned num_baryon = 0;

  // vector of selected clusters
  std::vector<CluVecIt>  sel;

  unsigned number_of_tries = _stepFactor*cv.size()*cv.size();
  if (number_of_tries<1) number_of_tries=1;

  long (*p_irnd)(long) = UseRandom::irnd;
  for (unsigned reconnections_tries = 0; reconnections_tries < number_of_tries; reconnections_tries++) {
    num_meson = 0;
    num_baryon = 0;

    // flag if we are able to find a suitable combinations of clusters
    bool _found = false;

    // Shuffle list of clusters to avoid systematic bias in cluster selection
    random_shuffle(newcv.begin(), newcv.end(), p_irnd);

    // loop over clustervector to find CR candidates
    for (CluVecIt cit = newcv.begin(); cit != newcv.end(); ++cit) {

      // skip the clusters to be deleted later from 3->2 cluster CR
      if (find(deleted.begin(), deleted.end(), *cit) != deleted.end()) continue;

      // avoid clusters already containing diuarks
      if (hasDiquark(*cit)) continue;

      // add to selection
      sel.push_back(cit);

      if (_clustersFarApart(sel)) {
        // reject far appart CR
        sel.pop_back();
        continue;
      }

      bool isMeson=((*cit)->numComponents() == 2);

      if ( isMeson && (num_meson ==0|| num_meson==1) && num_baryon ==0) {
        num_meson++;
        /**
         * now we habe either 1 or 2 mesonic clusters and have to continue
         */
        continue;
      } else if ( isMeson && (num_baryon == 1 || num_meson ==2)) {
        num_meson++;
        _found = true;
        /**
         * we have either 3 mesonic or 1 mesonic and 1 baryonic cluster
         * and try to colour reconnect
         */
        break;
      } else if (num_baryon ==0 && num_meson==0) {
        num_baryon++;
        /**
         * now we have 1 baryonic cluster and have to continue
         */
        continue;
      } else if (num_meson == 2) {
        /**
         * we already have 2 mesonic clusters and dont want a baryonic one
         * since this is an invalid selection
         */
        // remove previously added cluster
        sel.pop_back();
        continue;
      } else {
        num_baryon++;
        _found = true;
        /**
         * now we have either 2 baryonic clusters or 1 mesonic and 1 baryonic cluster
         * and try to colour reconnect
         */
        break;
      }
    }

    // added for more efficent rejection if some reco probabilities are 0
    if ( _found ) {

      // reject MBtoMB candidates if _precoMB_MB=0
      if ( _precoMB_MB == 0 && (num_baryon == 1 && num_meson == 1) ) {
        _found=false;
      }

      // reject BbarBto3M candidates if _precoBbarB_3M=0
      if ( _precoBbarB_3M== 0 && num_baryon == 2 ) {
        bool isBbarBto3Mcandidate=(
                                    (*sel[0])->particle(0)->hasColour() && (*sel[1])->particle(0)->hasColour(true) )
                                  || ( (*sel[0])->particle(0)->hasColour(true) && (*sel[1])->particle(0)->hasColour() );

        if ( isBbarBto3Mcandidate) _found=false;
      }

      // reject 2Bto2B candidates if _preco2B_2B=0
      if ( _preco2B_2B == 0 && num_baryon == 2 ) {
        bool is2Bto2Bcandidate=(
                                 (*sel[0])->particle(0)->hasColour() && (*sel[1])->particle(0)->hasColour() )
                               || ( (*sel[0])->particle(0)->hasColour(true) && (*sel[1])->particle(0)->hasColour(true) );

        if ( is2Bto2Bcandidate ) _found=false;
      }
    }
    // were we able to find a combination?
    if (_found==false) {
      // clear the selection if we did not find a valid set of  clusters
      sel.erase(sel.begin(), sel.end());
      continue;
    }
    assert(sel.size()<4);
    assert(sel.size()>1);

    string kind_of_reco = "";
    int reco_info[3];

    // find best CR option for the selection
    _findbestreconnectionoption(sel, num_baryon, kind_of_reco, reco_info);

    if (kind_of_reco == "") {
      // no reconnection was found
      sel.erase(sel.begin(), sel.end());
      continue;
    } else if (kind_of_reco == "3Mto3M" && UseRandom::rnd() < _preco3M_3M) {
      // 3Mto3M colour reconnection
      const auto & reconnected = _reconnect3Mto3M(*sel[0], *sel[1], *sel[2],
                                          reco_info);
      (*sel[0]) = std::get<0>(reconnected);
      (*sel[1]) = std::get<1>(reconnected);
      (*sel[2]) = std::get<2>(reconnected);
    } else if (kind_of_reco=="2Bto3M" && UseRandom::rnd() < _precoBbarB_3M) {
      // antibaryonic and baryonic to 3 mesonic reconnecion
      const auto & reconnected = _reconnectBbarBto3M(*sel[0], *sel[1],
                                             reco_info[0], reco_info[1], reco_info[2]);
      (*sel[0]) = std::get<0>(reconnected);
      (*sel[1]) = std::get<1>(reconnected);
      newcv.push_back(std::get<2>(reconnected));
    } else if (kind_of_reco=="3Mto2B"
				&& _canMakeBaryonicCluster(*sel[0], *sel[1], *sel[2])
				&& UseRandom::rnd() < _preco3M_BBbar) {
      // 3 mesonic to antibaryonic and baryonic reconnection
      ClusterPtr b1, b2;
      _makeBaryonicClusters(*sel[0], *sel[1], *sel[2], b1, b2);
      (*sel[0]) = b1;
      (*sel[1]) = b2;
      deleted.push_back(*sel[2]);
    } else if (kind_of_reco=="2Bto2B" && UseRandom::rnd() < _preco2B_2B) {
      // 2 (anti)baryonic to 2 (anti)baryonic reconnection
      const auto & reconnected = _reconnect2Bto2B(*sel[0], *sel[1],
                                          reco_info[0], reco_info[1]);
      (*sel[0]) = reconnected.first;
      (*sel[1]) = reconnected.second;
    } else if (kind_of_reco=="MBtoMB" && UseRandom::rnd() < _precoMB_MB) {
      // (anti)baryonic and mesonic to (anti)baryonic and mesonic reconnection
      const auto & reconnected = _reconnectMBtoMB(*sel[0], *sel[1],
                                          reco_info[0]);
      (*sel[0]) = reconnected.first;
      (*sel[1]) = reconnected.second;
    }
    // erase the sel-vector
    sel.erase(sel.begin(), sel.end());
  }

  // write to clustervector new CR'd clusters and deleting
  // all deleted clusters
  ClusterVector clustervector;
  clustervector.reserve(newcv.size());
  for (const auto & cluster : newcv)
    if (find(deleted.begin(), deleted.end(), cluster) == deleted.end())
      clustervector.push_back(cluster);

  swap(cv, clustervector);
}

namespace {
bool doublyHeavy(CluVecIt c1, CluVecIt c2) {
  if(abs((**c1).    colParticle()->id())>=4 &&
     abs((**c1).    colParticle()->id())>=4) return true;
  if(abs((**c1).antiColParticle()->id())>=4 &&
     abs((**c2).antiColParticle()->id())>=4) return true;
  return true;
}
}

void ColourReconnector::_findbestreconnectionoption(
    std::vector<CluVecIt> & cls, const unsigned & baryonic,
    string & kind_of_reco, int (&reco_info)[3]) const {
  double min_displacement;
  if (baryonic == 0) {
    // case with 3 mesonic clusters
    assert(cls.size() == 3);

    // calculate the initial displacement sum
    min_displacement  = _mesonToBaryonFactor * _displacement((*cls[0])->particle(0), (*cls[0])->particle(1));
    min_displacement += _mesonToBaryonFactor * _displacement((*cls[1])->particle(0), (*cls[1])->particle(1));
    min_displacement += _mesonToBaryonFactor * _displacement((*cls[2])->particle(0), (*cls[2])->particle(1));

    // find best CR reco_info and kind_of_reco
    _3MtoXreconnectionfinder(cls,
                             reco_info[0], reco_info[1], reco_info[2], min_displacement, kind_of_reco);
    /**
     * kind_of_reco either "3Mto3M" or "3Mto2B" (or "" if no better configuration is found)
     * case 3Mto3M:    the coloured particle of the i-th cluster forms a new cluster with the
     * 				   antiparticle of the reco_info[i]-th cluster
     * case 3MtoBbarB: all 3 (anti)coloured particle form a new (anti)baryonic cluster
     */
  } else if (baryonic == 1) {
    // case 1 baryonic and 1 mesonic cluster
    assert(cls.size() == 2);

    // make mesonic cluster always the cls[0]
    if ((*cls[0])->numComponents() == 3) {
      ClusterPtr zw = *cls[0];
      *cls[0] = *cls[1];
      *cls[1] = zw;
    }

    // calculate the initial displacement sum
    min_displacement  = _mesonToBaryonFactor *_displacement((*cls[0])->particle(0), (*cls[0])->particle(1));
    min_displacement += _displacementBaryonic((*cls[1])->particle(0), (*cls[1])->particle(1), (*cls[1])->particle(2));

    // find best CR reco_info and kind_of_reco
    _BMtoBMreconnectionfinder(*cls[0], *cls[1],
                              reco_info[0], min_displacement, kind_of_reco);
    /**
     * reco_info[0] is the index of the (anti)quarks of the baryonic cluster cls[1], which should
     * be swapped with the (anti)quarks of the mesonic cluster cls[0]
     */

  } else {
    assert(baryonic==2);
    assert(cls.size()==2);

    // calculate the initial displacement sum
    min_displacement  = _displacementBaryonic((*cls[0])->particle(0), (*cls[0])->particle(1), (*cls[0])->particle(2));
    min_displacement += _displacementBaryonic((*cls[1])->particle(0), (*cls[1])->particle(1), (*cls[1])->particle(2));

    // case 2 (anti)baryonic clusters to 2 other (anti)baryonic clusters
    if (      ( (*cls[0])->particle(0)->hasColour()     && (*cls[1])->particle(0)->hasColour()     )
              || ( (*cls[0])->particle(0)->hasColour(true) && (*cls[1])->particle(0)->hasColour(true) ) ) {
      // find best CR reco_info and kind_of_reco
      _2Bto2BreconnectionFinder(*cls[0], *cls[1],
                                reco_info[0], reco_info[1], min_displacement, kind_of_reco);
      /**
       * swap the reco_info[0]-th particle of the first cluster in the vector with the
       * reco_info[1]-th particle of the second cluster
       */
    } else {
      // case 1 baryonic and 1 antibaryonic cluster to 3 mesonic clusters

      // find best CR reco_info and kind_of_reco
      _BbarBto3MreconnectionFinder(*cls[0], *cls[1],
                                   reco_info[0], reco_info[1], reco_info[2], min_displacement, kind_of_reco);
      /**
       * the i-th particle of the first cluster form a new mesonic cluster with the
       * reco_info[i]-th particle of the second cluster
       */
    }
  }
  return;
}


void ColourReconnector::_findPartnerBaryonicDiquarkCluster(
  const CluVecIt & cl, ClusterVector & cv,
  unsigned & typeOfReconnection,
  const ClusterVector& deleted,
  CluVecIt &candidate1,
  CluVecIt &candidate2 ) const {

  typeOfReconnection = 0; // no Reconnection found
  using Constants::pi;
  using Constants::twopi;

  candidate1 = cl;
  candidate2 = cl;

  bool candIsOctet1 = false;
  bool candIsOctet2 = false;

  bool candIsQQ1 = false;
  bool candIsQQ2 = false;

  bool foundCR = false;

  double maxrap1 = 0.0;
  double maxrap2 = 0.0;

  double minrap1 = 0.0;
  double minrap2 = 0.0;
  
  double maxsum1 = 0.0;
  double maxsum2 = 0.0;

  double NegativeRapidtyThreshold = 0.0;
  double PositiveRapidtyThreshold = 0.0;

  // Get n and nbar direction for rapidity calculation
	const Boost pCluBoost = (*cl)->momentum().boostVector();
	Lorentz5Momentum p1Vec = (*cl)->antiColParticle()->momentum();
	p1Vec.boost(-pCluBoost);

	// n and nbar in COM frame of cluster and boost back
  LorentzVector<double> p1_n   ( p1Vec.vect().unit(),1);
  LorentzVector<double> p1_nbar(-p1Vec.vect().unit(),1);
	p1_n.boost(pCluBoost);
	p1_nbar.boost(pCluBoost);

  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {
    //avoid looping over clusters containing diquarks
    if ( (*cit)->numComponents()>=3 ) continue;

    //skip the cluster to be deleted later 3->2 cluster
    if ( find(deleted.begin(), deleted.end(), *cit) != deleted.end() )
      continue;

    if ( hasDiquark(*cit) ) continue;
    if ( (*cl)->isBeamCluster() && (*cit)->isBeamCluster() )
      continue;
    if ( cit==cl ) continue;

	// veto if Clusters further apart than _maxDistance
	if (_localCR && ((**cl).vertex().vect()-(**cit).vertex().vect()).mag() > _maxDistance) continue;
	// veto if Clusters have negative spacetime difference
	if (_causalCR && ((**cl).vertex()-(**cit).vertex()).m() < ZERO) continue;

	bool octetNormalCR = 
		    (_isColour8( (*cl)->colParticle(), (*cit)->antiColParticle() )
			||
			_isColour8( (*cit)->colParticle(), (*cl)->antiColParticle() ) );

    // boost constituents of cit into RF of cl
    Lorentz5Momentum p2col = (*cit)->colParticle()->momentum();
    Lorentz5Momentum p2anticol = (*cit)->antiColParticle()->momentum();

		// Here we compute the relative rapidity with respect to p1_n p1_nbar direction
    const double rapq = calculateRelativeRapidity(p1_n, p1_nbar, p2col);
    const double rapqbar = calculateRelativeRapidity(p1_n, p1_nbar, p2anticol);

	// configuration for Mesonic CR
	if ( rapq > 0.0 && rapqbar < 0.0
			&& rapq    > PositiveRapidtyThreshold
			&& rapqbar < NegativeRapidtyThreshold) {
		//sum of rapidities of quarks
		const double sumQQbar = abs(rapq) + abs(rapqbar);
		if ( sumQQbar > maxsum2 ) {
			if ( sumQQbar > maxsum1 ) {
				double factor = candIsQQ1 ? _mesonToBaryonFactor:1.0;
				maxsum2 = (factor*maxsum1) > sumQQbar ? sumQQbar:(factor*maxsum1);
				candidate2 = candidate1;
				candIsQQ2 = candIsQQ1;
				candIsOctet2 = candIsOctet1;
				maxrap2 = maxrap1;
				minrap2 = minrap1;

				maxsum1 = sumQQbar;
				candidate1 = cit;
				candIsQQ1 = false;
				candIsOctet1 = octetNormalCR;
				maxrap1 = rapq;
				minrap1 = rapqbar;
			} else {
				maxsum2 = sumQQbar;
				candidate2 = cit;
				candIsQQ2 = false;
				candIsOctet2 = octetNormalCR;
				maxrap2 = rapq;
				minrap2 = rapqbar;
			}
			// choose the less stringent threshold for further iterations
			PositiveRapidtyThreshold = maxrap1 > maxrap2 ? maxrap2:maxrap1;
			NegativeRapidtyThreshold = minrap1 < minrap2 ? minrap2:minrap1;
			foundCR=true;
		}
	}
	assert(PositiveRapidtyThreshold<=maxrap1);
	assert(PositiveRapidtyThreshold<=maxrap2);
	assert(NegativeRapidtyThreshold>=minrap1);
	assert(NegativeRapidtyThreshold>=minrap2);
	assert(maxsum1>=maxsum2);
	if ( rapq < 0.0 && rapqbar > 0.0
			&& rapqbar > PositiveRapidtyThreshold/_mesonToBaryonFactor
			&& rapq    < NegativeRapidtyThreshold/_mesonToBaryonFactor ) {
		//sum of rapidities of quarks
		const double sumQQ = abs(rapq) + abs(rapqbar);
		if ( sumQQ > maxsum2/_mesonToBaryonFactor ) {
			if ( sumQQ > maxsum1 ) {
				double factor = candIsQQ1 ? _mesonToBaryonFactor:1.0;
				maxsum2 = (factor*maxsum1) > sumQQ ? sumQQ:(factor*maxsum1);
				candidate2 = candidate1;
				candIsQQ2 = candIsQQ1;
				candIsOctet2 = candIsOctet1;
				maxrap2 = maxrap1;
				minrap2 = minrap1;

				maxsum1 = sumQQ;
				candidate1 = cit;
				candIsQQ1 = true;
				candIsOctet1 = octetNormalCR;
				maxrap1 = rapqbar;
				minrap1 = rapq;
			} else {
				maxsum2 = (_mesonToBaryonFactor*maxsum1) > sumQQ ? sumQQ:(_mesonToBaryonFactor*maxsum1);
				candidate2 = cit;
				candIsQQ2 = true;
				candIsOctet2 = octetNormalCR;
				maxrap2 = rapqbar;
				minrap2 = rapq;
			}
			// choose the less stringent threshold for further iterations
			PositiveRapidtyThreshold = maxrap1 > maxrap2 ? maxrap2:maxrap1;
			NegativeRapidtyThreshold = minrap1 < minrap2 ? minrap2:minrap1;
			foundCR=true;
		}
    }
	assert(PositiveRapidtyThreshold<=maxrap1);
	assert(PositiveRapidtyThreshold<=maxrap2);
	assert(NegativeRapidtyThreshold>=minrap1);
	assert(NegativeRapidtyThreshold>=minrap2);
	assert(maxsum1>=maxsum2);

  }
  // determine the type
  if (!candIsQQ1) {
		// Mesonic CR
	  if (candIsOctet1) {
			if (!candIsQQ2 && !candIsOctet2) {
				if (candidate2!=cl){
					swap(candidate2,candidate1);
					typeOfReconnection = 1;
				}
				else typeOfReconnection = 0;
			}
			else
				typeOfReconnection = 0;
		}
	  else
		  typeOfReconnection = 1;
  }
  else if (candIsQQ1)
  {
	  if (candIsQQ2 && _canMakeBaryonicCluster(*cl, *candidate1, *candidate2))
			// Baryonic CR
		  typeOfReconnection = 2;
	  else if (_canMakeDiquarkCluster(*cl, *candidate1))
			// Diquark CR
		  typeOfReconnection = 3;
		else
			// No CR
		  typeOfReconnection = 0;
  }
  if (!foundCR) typeOfReconnection = 0;
  // veto reconnection if cannot make a Diquark Cluster
  bool failDCR=false;
  if (typeOfReconnection == 3) {
	  if (_diquarkCR && !_canMakeDiquarkCluster(*cl, *candidate1)) {
			  // No CR
			  typeOfReconnection = 0;
			  failDCR=true;
	  }
  }
	if (typeOfReconnection == 2 && !_canMakeBaryonicCluster(*cl,*candidate1,*candidate2)) {
		std::cout << "WARNING Should never happen!!!!!" << std::endl;
		if (_canMakeDiquarkCluster(*cl,*candidate1)){
			// if we cannot make baryonic CR try diquark CR
			typeOfReconnection = 3;
		}
		else {
			// reject CR if we cannot make neither baryonic nor diquark CR
			typeOfReconnection = 0;
		}
	}
  if (_debug) {
	  	std::ofstream outTypes("WriteOut/TypesOfDCR.dat", std::ios::app);
		outTypes << (failDCR ? 4:typeOfReconnection) << "\n";
		outTypes.close();
		  switch (typeOfReconnection)
		  {
			  // Mesonic CR
			  case 1:
				  {
				  std::ofstream outMCR("WriteOut/MCR.dat", std::ios::app);
				  outMCR << minrap1 << "\t"
						 << maxrap1 << "\t"
						 << minrap2 << "\t"
						 << maxrap2 << "\t"
						 <<"\n";
				  outMCR.close();
				  break;
				  }
				  // Baryonic CR
			  case 2:
				  {
				  std::ofstream outBCR("WriteOut/BCR.dat", std::ios::app);
				  outBCR << minrap1 << "\t"
						 << maxrap1 << "\t"
						 << minrap2 << "\t"
						 << maxrap2 << "\t"
						 <<"\n";
				  outBCR.close();
				  break;
				  }
				  // Diquark CR
			  case 3:
				  {
				  std::ofstream outDCR("WriteOut/DCR.dat", std::ios::app);
				  outDCR << minrap1 << "\t"
						 << maxrap1 << "\t"
						 << minrap2 << "\t"
						 << maxrap2 << "\t"
						 <<"\n";
				  outDCR.close();
				  break;
				  }
				  // No CR found
			  case 0:
				  {
				  std::ofstream outNoCR("WriteOut/NoCR.dat", std::ios::app);
				  outNoCR<< minrap1 << "\t"
						 << maxrap1 << "\t"
						 << minrap2 << "\t"
						 << maxrap2 << "\t"
						 <<"\n";
				  outNoCR.close();
				  break;
				  }
			  default:
				  assert(false);
		  }
  }
}


CluVecIt ColourReconnector::_findPartnerBaryonic(
  const CluVecIt & cl, ClusterVector & cv,
  bool & baryonicCand,
  const ClusterVector& deleted,
  CluVecIt &baryonic1,
  CluVecIt &baryonic2 ) const {

  using Constants::pi;
  using Constants::twopi;

  // Returns a candidate for possible reconnection
  CluVecIt candidate = cl;

  bool bcand = false;

  double maxrap = 0.0;
  double minrap = 0.0;
  double maxrapNormal = 0.0;
  double minrapNormal = 0.0;
  double maxsumnormal = 0.0;

  double maxsum = 0.0;
  double secondsum = 0.0;

  // Get n and nbar direction for rapidity calculation
	const Boost pCluBoost = (*cl)->momentum().boostVector();
	Lorentz5Momentum p1Vec = (*cl)->antiColParticle()->momentum();
	p1Vec.boost(-pCluBoost);

	// n and nbar in COM frame of cluster and boost back
  LorentzVector<double> p1_n   ( p1Vec.vect().unit(),1);
  LorentzVector<double> p1_nbar(-p1Vec.vect().unit(),1);
	p1_n.boost(pCluBoost);
	p1_nbar.boost(pCluBoost);

  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {
    //avoid looping over clusters containing diquarks
    if ( hasDiquark(*cit) ) continue;
    if ( (*cit)->numComponents()>=3 ) continue;
    if ( cit==cl ) continue;

    //skip the cluster to be deleted later 3->2 cluster
    if ( find(deleted.begin(), deleted.end(), *cit) != deleted.end() )
      continue;

    if ( (*cl)->isBeamCluster() && (*cit)->isBeamCluster() )
      continue;

    // veto if Clusters further apart than _maxDistance
    if (_localCR && ((**cl).vertex().vect()-(**cit).vertex().vect()).mag() > _maxDistance) continue;
    // veto if Clusters have negative spacetime difference
    if (_causalCR && ((**cl).vertex()-(**cit).vertex()).m() < ZERO) continue;

    const bool Colour8 =
      _isColour8( (*cl)->colParticle(), (*cit)->antiColParticle() )
      ||
      _isColour8( (*cit)->colParticle(), (*cl)->antiColParticle() ) ;


    // boost constituents of cit into RF of cl
    Lorentz5Momentum p2col = (*cit)->colParticle()->momentum();
    Lorentz5Momentum p2anticol = (*cit)->antiColParticle()->momentum();

    // Here we compute the relative rapidity with respect to p1_n p1_nbar direction
    const double rapq = calculateRelativeRapidity(p1_n, p1_nbar, p2col);
    const double rapqbar = calculateRelativeRapidity(p1_n, p1_nbar, p2anticol);

    // configuration for normal CR
    if ( !Colour8
        && rapq > 0.0 && rapqbar < 0.0
        && rapq > maxrap
        && rapqbar < minrap ) {
      maxrap = rapq;
      minrap = rapqbar;
      //sum of rapidities of quarks
      const double normalsum = abs(rapq) + abs(rapqbar);
      if ( normalsum > maxsumnormal ) {
        maxsumnormal = normalsum;
        maxrapNormal = rapq;
        minrapNormal = rapqbar;
        bcand = false;
        candidate = cit;
      }

    }
    if ( rapq < 0.0 && rapqbar >0.0
        && rapqbar > maxrapNormal
        && rapq < minrapNormal ) {
      maxrap = rapqbar;
      minrap = rapq;
      const double sumrap = abs(rapqbar) + abs(rapq);
      // first candidate gets here. If second baryonic candidate has higher Ysum than the first
      // one, the second candidate becomes the first one and the first the second.
      if (sumrap > maxsum) {
        if (maxsum != 0) {
	  if(_doublyHeavyBaryon || (!doublyHeavy(candidate,cit)&&
				    !doublyHeavy(baryonic1,cit))) {
	    baryonic2 = baryonic1;
	    baryonic1 = cit;
	    bcand = true;
	  }
	  else if(!doublyHeavy(candidate,cit)) {
	    baryonic1 = cit;
	  }
        }
	else {
	  if(_doublyHeavyBaryon || !doublyHeavy(candidate,cit))
	    baryonic1 = cit;
        }
        maxsum = sumrap;
      }
      else {
        if (sumrap > secondsum && sumrap != maxsum) {
	  if(_doublyHeavyBaryon || (!doublyHeavy(candidate,cit)&&
				    !doublyHeavy(baryonic1,cit))) {
	    secondsum = sumrap;
	    bcand = true;
	    baryonic2 = cit;
	  }
        }
      }
    }
  }

  if (bcand == true) {
    baryonicCand = true;
  }

  if (_debug) {
	  	std::ofstream outTypes("WriteOut/TypesOfBCR.dat", std::ios::app);
		outTypes << (baryonicCand ? 2:(candidate==cl ? 0:1)) << "\n";
		outTypes.close();
  }
  return candidate;
}

CluVecIt ColourReconnector::_findRecoPartnerPlainDynamic(const CluVecIt & cl,
    ClusterVector & cv, const ClusterVector & deleted, bool diquarkCR) const {

  CluVecIt candidate = cl;
	double pNoRecoMin=1.0;
	double pNoReco,preco,precoDiquark;

  // Get n and nbar direction for rapidity calculation
	const Boost pCluBoost = (*cl)->momentum().boostVector();
	Lorentz5Momentum p1Vec = (*cl)->antiColParticle()->momentum();
	p1Vec.boost(-pCluBoost);

	// n and nbar in COM frame of cluster and boost back
  LorentzVector<double> p1_n   ( p1Vec.vect().unit(),1);
  LorentzVector<double> p1_nbar(-p1Vec.vect().unit(),1);
	p1_n.boost(pCluBoost);
	p1_nbar.boost(pCluBoost);

  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {

		// Skip deleted clusters
		if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
			continue;
		// skip diquark clusters
    if ((*cit)->numComponents()>2 || (_diquarkCR && hasDiquark(*cit))) continue;

    // don't even look at original cluster
    if (cit==cl) continue;

    // don't allow colour octet clusters
    bool Colour8 = ( _isColour8( (*cl)->colParticle(),
                     (*cit)->antiColParticle() )  ||
         _isColour8( (*cit)->colParticle(),
                     (*cl)->antiColParticle() ) );
    // stop it putting beam remnants together
    if ((*cl)->isBeamCluster() && (*cit)->isBeamCluster()) continue;

		// veto if Clusters further apart than _maxDistance
		if (_localCR && ((**cl).vertex().vect()-(**cit).vertex().vect()).mag() > _maxDistance) continue;
		// veto if Clusters have negative spacetime difference
		if (_causalCR && ((**cl).vertex()-(**cit).vertex()).m() < ZERO) continue;

    // boost constituents of cit into RF of cl
    Lorentz5Momentum p2col = (*cit)->colParticle()->momentum();
    Lorentz5Momentum p2anticol = (*cit)->antiColParticle()->momentum();

		// Here we compute the relative rapidity with respect to p1_n, p1_nbar direction
    const double rapq = calculateRelativeRapidity(p1_n, p1_nbar, p2col);
    const double rapqbar = calculateRelativeRapidity(p1_n, p1_nbar, p2anticol);

    // configuration for normal CR
		if ( !Colour8
				&& rapq > 0.0 && rapqbar < 0.0) {
		}
    // configuration for diquark CR
		else if ( _diquarkCR && rapq < 0.0 && rapqbar > 0.0) {

		}
		else
			continue;

		// Get dynamic CR probabilities and try to minimize the non-reco probability
		std::tie( pNoReco, preco, precoDiquark) = _dynamicRecoProbabilitiesCF2(*cl,*cit,diquarkCR);
		if (pNoReco<pNoRecoMin) {
			candidate = cit;
			pNoRecoMin=pNoReco;
		}
  }

	if (_debug) {
		ofstream out("WriteOut/pNoReco.dat", std::ios::app);
		out << pNoRecoMin << "\t" << cv.size() << "\n";
		out.close();
	}
  return candidate;
}


CluVecIt ColourReconnector::_findRecoPartnerPlain(const CluVecIt & cl,
    ClusterVector & cv, const ClusterVector & deleted) const {

  CluVecIt candidate = cl;
  Energy minMass = 1*TeV;
  for (CluVecIt cit=cv.begin(); cit != cv.end(); ++cit) {

		// Skip deleted clusters
		if (find(deleted.begin(), deleted.end(), *cit) != deleted.end())
			continue;
		// skip diquark clusters
    if ((*cit)->numComponents()>2 || (_diquarkCR && hasDiquark(*cit))) continue;


    // don't even look at original cluster
    if (cit==cl) continue;

    // don't allow colour octet clusters
    if ( _isColour8( (*cl)->colParticle(),
                     (*cit)->antiColParticle() )  ||
         _isColour8( (*cit)->colParticle(),
                     (*cl)->antiColParticle() ) ) {
      continue;
    }

    // stop it putting beam remnants together
    if ((*cl)->isBeamCluster() && (*cit)->isBeamCluster()) continue;

	// veto if Clusters further apart than _maxDistance
	if (_localCR && ((**cl).vertex().vect()-(**cit).vertex().vect()).mag() > _maxDistance) continue;
	// veto if Clusters have negative spacetime difference
	if (_causalCR && ((**cl).vertex()-(**cit).vertex()).m() < ZERO) continue;

    // momenta of the old clusters
    Lorentz5Momentum p1 = (*cl)->colParticle()->momentum() +
                          (*cl)->antiColParticle()->momentum();
    Lorentz5Momentum p2 = (*cit)->colParticle()->momentum() +
                          (*cit)->antiColParticle()->momentum();

    // momenta of the new clusters
    Lorentz5Momentum p3 = (*cl)->colParticle()->momentum() +
                          (*cit)->antiColParticle()->momentum();
    Lorentz5Momentum p4 = (*cit)->colParticle()->momentum() +
                          (*cl)->antiColParticle()->momentum();

		Energy oldMass = abs( p1.m() ) + abs( p2.m() );
		Energy newMass = abs( p3.m() ) + abs( p4.m() );


		if ( newMass < oldMass && newMass < minMass ) {
			minMass = newMass;
			candidate = cit;
		}
  }

  return candidate;
}

// forms two baryonic clusters from three clusters
void ColourReconnector::_makeBaryonicClusters(
  ClusterPtr &c1, ClusterPtr &c2,
  ClusterPtr &c3,
  ClusterPtr &newcluster1,
  ClusterPtr &newcluster2) const {

  //make sure they all have 2 components
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  assert(c3->numComponents()==2);
  // check if doubly heavy
  if( (int(abs(c1->    colParticle()->id())>=4) + int(abs(c2->    colParticle()->id())>=4)+ int(abs(c3->    colParticle()->id())>=4) >=2) ||
      (int(abs(c1->antiColParticle()->id())>=4) + int(abs(c2->antiColParticle()->id())>=4)+ int(abs(c3->antiColParticle()->id())>=4) >=2) )
    generator()->log() << "Colour recon making doubly heavy baryon "
		       << *c1->    colParticle() << "\n" 
		       << *c2->    colParticle() << "\n" 
		       << *c3->    colParticle() << "\n" 
		       << *c1->antiColParticle() << "\n" 
		       << *c2->antiColParticle() << "\n" 
		       << *c3->antiColParticle() << "\n"; 
  //abandon children
  c1->colParticle()->abandonChild(c1);
  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);
  c2->antiColParticle()->abandonChild(c2);
  c3->colParticle()->abandonChild(c3);
  c3->antiColParticle()->abandonChild(c3);

  newcluster1 = new_ptr(Cluster(c1->colParticle(),c2->colParticle(), c3->colParticle()));
  c1->colParticle()->addChild(newcluster1);
  c2->colParticle()->addChild(newcluster1);
  c3->colParticle()->addChild(newcluster1);
  newcluster1->setVertex(LorentzPoint());

  newcluster2 = new_ptr(Cluster(c1->antiColParticle(), c2->antiColParticle(),
                                c3->antiColParticle()));
  c1->antiColParticle()->addChild(newcluster2);
  c2->antiColParticle()->addChild(newcluster2);
  c3->antiColParticle()->addChild(newcluster2);
  newcluster2->setVertex(LorentzPoint());
}
bool ColourReconnector::_canMakeDiquarkCluster(tcPPtr pCol1, tcPPtr pCol2,tcPPtr pAntiCol1, tcPPtr pAntiCol2) const{
	double a;
	return _canMakeDiquarkCluster( pCol1,  pCol2, pAntiCol1,  pAntiCol2,a);
}
bool ColourReconnector::_canMakeDiquarkCluster(tcPPtr pCol1, tcPPtr pCol2,tcPPtr pAntiCol1, tcPPtr pAntiCol2, double & PhaseSpace) const{
	tcPDPtr dataDiquark    =  _hadronSpectrum->makeDiquark(pCol1->dataPtr(), pCol2->dataPtr());
	tcPDPtr dataDiquarkBar =  _hadronSpectrum->makeDiquark(pAntiCol1->dataPtr(), pAntiCol2->dataPtr());
	if (!dataDiquark){
		throw Exception() << "Could not make a diquark from"
			<< pCol1->dataPtr()->PDGName() << " and "
			<< pCol2->dataPtr()->PDGName()
			<< " in ColourReconnector::_canMakeDiquarkCluster()"
			<< Exception::eventerror;
	}
	if (!dataDiquarkBar){
		throw Exception() << "Could not make an anti-diquark from"
			<< pAntiCol1->dataPtr()->PDGName() << " and "
			<< pAntiCol2->dataPtr()->PDGName()
			<< " in ColourReconnector::_canMakeDiquarkCluster()"
			<< Exception::eventerror;
	}
	int diqTreatment = _clusterFinder->diquarkOnShell();
	Lorentz5Momentum Ptot=pCol1->momentum() + pCol2->momentum() + pAntiCol1->momentum() + pAntiCol2->momentum();
	Energy Mass=Ptot.m();
	Energy minMassCD = _hadronSpectrum->massLightestHadronPair(dataDiquark,dataDiquarkBar);
	// We need to guarantee that a diquark cluster can decay to the two lightest hadrons
	if ( Mass<=minMassCD ) {
		// DiQuarkOnShell = Mixed (Default) -> diqTreatment==-1
		// DiQuarkOnShell = No -> diqTreatment==0
		return false;
	}
	
	if (diqTreatment == 1){
		// DiQuarkOnShell = Yes
		Energy minMassOnShell = dataDiquark->constituentMass() + dataDiquarkBar->constituentMass();

		// We need to guarantee that a diquark cluster can have its two future diquarks on shell
		if ( Mass<=minMassOnShell ) {
			return false;
		}
	}	

	if (_cutDiquarkClusterFormation>0) {
		double cut = _cutDiquarkClusterFormation;
		double valDiq = pCol1->momentum()*pCol2->momentum()/(pCol1->momentum().m()*pCol2->momentum().m())-1.0;
		double valAntiDiq = pAntiCol1->momentum()*pAntiCol2->momentum()/(pAntiCol1->momentum().m()*pAntiCol2->momentum().m())-1.0;
		if (valDiq>cut)
			return false;
		if (valAntiDiq>cut)
			return false;
	}

	if (_phaseSpaceDiquarkFission){
		// Tried to add a factor suppressing to low mass Diquark Clusters
		double factor;
		switch (_phaseSpaceDiquarkFission)
		{
			case 1:
				{
					const auto & BaryonPair=_hadronSpectrum->lightestHadronPair(dataDiquark,dataDiquarkBar);
					if (Mass-(BaryonPair.first->mass()+BaryonPair.second->mass())<=ZERO)
						return false;
					factor = 2.0*Kinematics::pstarTwoBodyDecay(Mass,BaryonPair.first->mass(),BaryonPair.second->mass())/Mass;	
					break;
				}
			case 2:
				{
					if (Mass-(dataDiquark->constituentMass()+dataDiquarkBar->constituentMass())<ZERO)
						return false;
					factor = 2.0*Kinematics::pstarTwoBodyDecay(Mass,dataDiquark->constituentMass(),dataDiquarkBar->constituentMass())/Mass;	
					break;
				}
			default:
				assert(false);
		}
		
		assert(factor>=0.0);
		assert(factor<=1.0);
		PhaseSpace = factor;
	}
	else {
		PhaseSpace = 1.0;
	}
	return true;
}
bool ColourReconnector::_canMakeBaryonicCluster(const ClusterPtr &c1, const ClusterPtr &c2, const ClusterPtr &c3) const {
    //make sure they all have 2 components
    assert(c1->numComponents()==2);
    assert(c2->numComponents()==2);
    assert(c3->numComponents()==2);
		vector<tcPPtr> Baryon = {c1->colParticle(), c2->colParticle(), c3->colParticle()};
		vector<tcPPtr> AntiBaryon = {c1->antiColParticle(), c2->antiColParticle(), c3->antiColParticle()};
		return (_canMakeBaryonicCluster(Baryon)
				&& _canMakeBaryonicCluster(AntiBaryon));
}

bool ColourReconnector::_canMakeBaryonicCluster(vector<tcPPtr> pCol) const {
	tcPDPtr dataDiquark12 = _hadronSpectrum->makeDiquark(pCol[0]->dataPtr(), pCol[1]->dataPtr());
	tcPDPtr dataDiquark13 = _hadronSpectrum->makeDiquark(pCol[0]->dataPtr(), pCol[2]->dataPtr());
	tcPDPtr dataDiquark23 = _hadronSpectrum->makeDiquark(pCol[1]->dataPtr(), pCol[2]->dataPtr());

	if (!dataDiquark12){
		throw Exception() << "Could not make a diquark from"
			<< pCol[0]->dataPtr()->PDGName() << " and "
			<< pCol[1]->dataPtr()->PDGName()
			<< " in ColourReconnector::_canMakeBaryonicCluster()"
			<< Exception::eventerror;
	}
	if (!dataDiquark13){
		throw Exception() << "Could not make a diquark from"
			<< pCol[0]->dataPtr()->PDGName() << " and "
			<< pCol[2]->dataPtr()->PDGName()
			<< " in ColourReconnector::_canMakeBaryonicCluster()"
			<< Exception::eventerror;
	}
	if (!dataDiquark23){
		throw Exception() << "Could not make a diquark from"
			<< pCol[1]->dataPtr()->PDGName() << " and "
			<< pCol[2]->dataPtr()->PDGName()
			<< " in ColourReconnector::_canMakeBaryonicCluster()"
			<< Exception::eventerror;
	}
	int diqTreatment = _clusterFinder->diquarkOnShell();
	Lorentz5Momentum Ptot=pCol[0]->momentum() + pCol[1]->momentum() + pCol[2]->momentum();
	Energy Mass=Ptot.m();
	Energy minMassCD = _hadronSpectrum->massLightestHadron(dataDiquark12, pCol[2]->dataPtr());
	if ( Mass<=minMassCD ) {
		// DiQuarkOnShell = Mixed (Default) -> diqTreatment==-1
		// DiQuarkOnShell = No -> diqTreatment==0
		return false;
	}
	
	if (diqTreatment == 1){
		// DiQuarkOnShell = Yes
		Energy minMassOnShell12 = dataDiquark12->constituentMass() + pCol[2]->dataPtr()->constituentMass();
		Energy minMassOnShell13 = dataDiquark13->constituentMass() + pCol[1]->dataPtr()->constituentMass();
		Energy minMassOnShell23 = dataDiquark23->constituentMass() + pCol[0]->dataPtr()->constituentMass();

		if ( Mass<=minMassOnShell12 || Mass<=minMassOnShell13 || Mass<=minMassOnShell23 ) {
			return false;
		}
	}
	return true;
}



bool ColourReconnector::_canMakeDiquarkCluster(const ClusterPtr &c1, const ClusterPtr &c2) const {
  double a;
  return _canMakeDiquarkCluster(c1,c2,a);
}

// forms a four-quark cluster
bool ColourReconnector::_canMakeDiquarkCluster(const ClusterPtr &c1, const ClusterPtr &c2, double & PhaseSpace) const {
  //make sure they all have 2 components
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  return _canMakeDiquarkCluster(c1->colParticle(),c2->colParticle(),c1->antiColParticle(),c2->antiColParticle(),PhaseSpace);
}
bool ColourReconnector::_makeDiquarkCluster(
                ClusterPtr &c1, ClusterPtr &c2,
                ClusterPtr &newcluster) const{
  if (!_canMakeDiquarkCluster(c1,c2)){
    std::cout << "Should never execute! check this earlier" << std::endl;
    return false;
  }

  //abandon children
  c1->colParticle()->abandonChild(c1);
  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);
  c2->antiColParticle()->abandonChild(c2);

  // Note: the convention that c1->newcluster.particle(0,2) and c2->newcluster.particle(1,3)
  newcluster = new_ptr(Cluster(c1->colParticle(),c2->colParticle(),
        c1->antiColParticle(), c2->antiColParticle()));
  c1->colParticle()->addChild(newcluster);
  c2->colParticle()->addChild(newcluster);
  c1->antiColParticle()->addChild(newcluster);
  c2->antiColParticle()->addChild(newcluster);
  newcluster->setVertex(LorentzPoint());
  return true;
}

pair <ClusterPtr,ClusterPtr> ColourReconnector::_splitDiquarkCluster(
    ClusterPtr &diquarkCluster, bool colourReconnect) const{
  ClusterPtr newcluster1,newcluster2;
  //abandon children
  diquarkCluster->particleB(0)->abandonChild(diquarkCluster);
  diquarkCluster->particleB(1)->abandonChild(diquarkCluster);
  diquarkCluster->particleB(2)->abandonChild(diquarkCluster);
  diquarkCluster->particleB(3)->abandonChild(diquarkCluster);

  // Note: the convention that c1->newcluster.particle(0,2) and c2->newcluster.particle(1,3)
  // Here we decide if we want to colour reconnect the original clusters (0,2) and (1,3) to
  // the clusters (0,3) and (1,2)
  if (colourReconnect) {
    // colour reconnect
    newcluster1 = new_ptr(Cluster(diquarkCluster->particleB(0), diquarkCluster->particleB(3)));
    newcluster2 = new_ptr(Cluster(diquarkCluster->particleB(1), diquarkCluster->particleB(2)));
    diquarkCluster->particleB(0)->addChild(newcluster1);
    diquarkCluster->particleB(1)->addChild(newcluster2);
    diquarkCluster->particleB(2)->addChild(newcluster2);
    diquarkCluster->particleB(3)->addChild(newcluster1);
  }
  else {
    // keep original colour connection
    newcluster1 = new_ptr(Cluster(diquarkCluster->particleB(0), diquarkCluster->particleB(2)));
    newcluster2 = new_ptr(Cluster(diquarkCluster->particleB(1), diquarkCluster->particleB(3)));
    diquarkCluster->particleB(0)->addChild(newcluster1);
    diquarkCluster->particleB(1)->addChild(newcluster2);
    diquarkCluster->particleB(2)->addChild(newcluster1);
    diquarkCluster->particleB(3)->addChild(newcluster2);
  }
  newcluster1->setVertex(LorentzPoint());
  newcluster2->setVertex(LorentzPoint());
  return pair <ClusterPtr, ClusterPtr> (newcluster1, newcluster2);
}

pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnect2Bto2B(ClusterPtr &c1, ClusterPtr &c2, const int s1, const int s2) const {

  // form the first new cluster

  // separate the quarks from their original cluster
  c1->particleB((s1+1)%3)->abandonChild(c1);
  c1->particleB((s1+2)%3)->abandonChild(c1);
  c2->particleB(s2)->abandonChild(c2);

  // now the new cluster
  ClusterPtr newCluster1 = new_ptr(Cluster(c1->particleB((s1+1)%3), c1->particleB((s1+2)%3), c2->particleB(s2)));

  c1->particleB((s1+1)%3)->addChild(newCluster1);
  c1->particleB((s1+2)%3)->addChild(newCluster1);
  c2->particleB(s2)->addChild(newCluster1);

  // set new vertex
  newCluster1->setVertex(LorentzPoint());

  // set beam remnants for new cluster
  if (c1->isBeamRemnant((s1+1)%3)) newCluster1->setBeamRemnant(0, true);
  if (c1->isBeamRemnant((s1+2)%3)) newCluster1->setBeamRemnant(1, true);
  if (c2->isBeamRemnant(s2)) newCluster1->setBeamRemnant(2, true);

  // for the second cluster same  procedure
  c2->particleB((s2+1)%3)->abandonChild(c2);
  c2->particleB((s2+2)%3)->abandonChild(c2);
  c1->particleB(s1)->abandonChild(c1);

  ClusterPtr newCluster2 = new_ptr(Cluster(c2->particleB((s2+1)%3), c2->particleB((s2+2)%3), c1->particleB(s1)));

  c2->particleB((s2+1)%3)->addChild(newCluster2);
  c2->particleB((s2+2)%3)->addChild(newCluster2);
  c1->particleB(s1)->addChild(newCluster2);

  newCluster2->setVertex(LorentzPoint());

  if (c2->isBeamRemnant((s2+1)%3)) newCluster2->setBeamRemnant(0, true);
  if (c2->isBeamRemnant((s2+2)%3)) newCluster2->setBeamRemnant(1, true);
  if (c1->isBeamRemnant(s1)) newCluster2->setBeamRemnant(2, true);

  return pair <ClusterPtr, ClusterPtr> (newCluster1, newCluster2);
}


std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr>
ColourReconnector::_reconnectBbarBto3M(ClusterPtr & c1, ClusterPtr & c2, const int s0, const int s1, const int s2) const {
  // make sure they all have 3 components
  assert(c1->numComponents()==3);
  assert(c2->numComponents()==3);

  // first Cluster
  c1->particleB(0)->abandonChild(c1);
  c2->particleB(s0)->abandonChild(c2);

  ClusterPtr newCluster1 = new_ptr(Cluster(c1->particleB(0), c2->particleB(s0)));

  c1->particleB(0)->addChild(newCluster1);
  c2->particleB(s0)->addChild(newCluster1);

  // set new vertex
  newCluster1->setVertex(0.5*(c1->particleB(0)->vertex() + c2->particleB(s0)->vertex()));

  // set beam remnants for new cluster
  if (c1->isBeamRemnant(0)) newCluster1->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

  // same for second cluster
  c1->particleB(1)->abandonChild(c1);
  c2->particleB(s1)->abandonChild(c2);

  ClusterPtr newCluster2 = new_ptr(Cluster(c1->particleB(1), c2->particleB(s1)));

  c1->particleB(1)->addChild(newCluster2);
  c2->particleB(s1)->addChild(newCluster2);

  newCluster2->setVertex(0.5*(c1->particleB(1)->vertex() + c2->particleB(s1)->vertex()));

  if (c1->isBeamRemnant(1)) newCluster2->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s1)) newCluster2->setBeamRemnant(1, true);

  // same for third cluster
  c1->particleB(2)->abandonChild(c1);
  c2->particleB(s2)->abandonChild(c2);

  ClusterPtr newCluster3 = new_ptr(Cluster(c1->particleB(2), c2->particleB(s2)));

  c1->particleB(2)->addChild(newCluster3);
  c2->particleB(s2)->addChild(newCluster3);

  newCluster3->setVertex(0.5*(c1->particleB(2)->vertex() + c2->particleB(s2)->vertex()));

  if (c1->isBeamRemnant(2)) newCluster3->setBeamRemnant(0, true);
  if (c2->isBeamRemnant(s2)) newCluster3->setBeamRemnant(1, true);

  return std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr> (newCluster1, newCluster2, newCluster3);
}

pair <ClusterPtr,ClusterPtr>
ColourReconnector::_reconnect(ClusterPtr &c1, ClusterPtr &c2) const {
	if (_becomesColour8Cluster(c1,c2)){
		std::cout << "Should never execute! _isColour8 in reconnect check this earlier" << std::endl;
	}

  // choose the other possibility to form two clusters from the given
  // constituents

  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  int c1_col(-1),c1_anti(-1),c2_col(-1),c2_anti(-1);
  for(unsigned int ix=0; ix<2; ++ix) {
    if     (c1->particle(ix)->hasColour(false)) c1_col  = ix;
    else if(c1->particle(ix)->hasColour(true )) c1_anti = ix;
    if     (c2->particle(ix)->hasColour(false)) c2_col  = ix;
    else if(c2->particle(ix)->hasColour(true )) c2_anti = ix;
  }
  assert(c1_col>=0&&c2_col>=0&&c1_anti>=0&&c2_anti>=0);

  c1->colParticle()->abandonChild(c1);
  c2->antiColParticle()->abandonChild(c2);

  ClusterPtr newCluster1
    = new_ptr( Cluster( c1->colParticle(), c2->antiColParticle() ) );

  newCluster1->setVertex(0.5*(c1->colParticle()->vertex() +
                              c2->antiColParticle()->vertex()));

  if(c1->isBeamRemnant(c1_col )) newCluster1->setBeamRemnant(0,true);
  if(c2->isBeamRemnant(c2_anti)) newCluster1->setBeamRemnant(1,true);

  c1->antiColParticle()->abandonChild(c1);
  c2->colParticle()->abandonChild(c2);

  ClusterPtr newCluster2
    = new_ptr( Cluster( c2->colParticle(), c1->antiColParticle() ) );


  c1->colParticle()->addChild(newCluster1);
  c2->antiColParticle()->addChild(newCluster1);
  c1->antiColParticle()->addChild(newCluster2);
  c2->colParticle()->addChild(newCluster2);

  newCluster2->setVertex(0.5*(c2->colParticle()->vertex() +
                              c1->antiColParticle()->vertex()));

  if(c2->isBeamRemnant(c2_col )) newCluster2->setBeamRemnant(0,true);
  if(c1->isBeamRemnant(c1_anti)) newCluster2->setBeamRemnant(1,true);

  return pair <ClusterPtr,ClusterPtr> (newCluster1, newCluster2);
}

std::tuple  <ClusterPtr, ClusterPtr> ColourReconnector::_reconnect3MtoMD(
    ClusterVector & cluvec, const int topology) const {
	assert(cluvec.size()==3);
	assert(topology>0 && topology<=33);
	int colIndexMCR = (topology/10)-1;
	int antiColIndexMCR = (topology-(colIndexMCR+1)*10)-1;
	assert(colIndexMCR<3 && colIndexMCR>=0);
	assert(antiColIndexMCR<3 && antiColIndexMCR>=0);

	ClusterPtr cMCR1 = cluvec[colIndexMCR];
	ClusterPtr cMCR2 = cluvec[antiColIndexMCR];


	if (_isColour8(cMCR1->colParticle(),cMCR2->antiColParticle())){
		std::cout << "Should never execute! _isColour8 in _reconnect3MtoMD check this earlier" << std::endl;
	}

	int colIdx[3]={-1,-1,-1};
	int anticolIdx[3]={-1,-1,-1};
	for (int i = 0; i < 3; i++) {
		assert(cluvec[i]->numComponents()==2);
		for (unsigned ip= 0; ip < 2; ip++){
			if (cluvec[i]->particle(ip)->hasColour(false)) colIdx[i]  = ip;
			if (cluvec[i]->particle(ip)->hasColour(true)) anticolIdx[i]  = ip;
		}
		assert(colIdx[i]>=0);
		assert(anticolIdx[i]>=0);
		// abbandon all children
		cluvec[i]->colParticle()->abandonChild(cluvec[i]);
		cluvec[i]->antiColParticle()->abandonChild(cluvec[i]);
	}
	// make the MCR cluster with indices colIndexMCR,antiColIndexMCR
	// form new mesonic cluster
	ClusterPtr newClusterMCR = new_ptr(Cluster(cluvec[colIndexMCR]->colParticle(), cluvec[antiColIndexMCR]->antiColParticle()));

	newClusterMCR->colParticle()->addChild(newClusterMCR);
	newClusterMCR->antiColParticle()->addChild(newClusterMCR);

	// set new vertex
	newClusterMCR->setVertex(0.5*(newClusterMCR->colParticle()->vertex() +
				newClusterMCR->antiColParticle()->vertex()));

	// set beam remnants for new cluster
	if (cluvec[colIndexMCR]->isBeamRemnant(colIdx[colIndexMCR])) newClusterMCR->setBeamRemnant(0, true);
	if (cluvec[antiColIndexMCR]->isBeamRemnant(anticolIdx[antiColIndexMCR])) newClusterMCR->setBeamRemnant(1, true);
	
	// make the DCR cluster with remaining indices
	int indexDCRcol1 = (colIndexMCR+1)%3;
	int indexDCRcol2 = (colIndexMCR+2)%3;
	int indexDCRacol1 = (antiColIndexMCR+1)%3;
	int indexDCRacol2 = (antiColIndexMCR+2)%3;
	assert(indexDCRcol1!=indexDCRcol2);
	assert(indexDCRcol1!=colIndexMCR);
	assert(indexDCRcol2!=colIndexMCR);
	assert(indexDCRacol1!=indexDCRacol2);
	assert(indexDCRacol1!=antiColIndexMCR);
	assert(indexDCRacol2!=antiColIndexMCR);
	ClusterPtr newClusterDCR = new_ptr(Cluster(
				cluvec[indexDCRcol1]->colParticle(),
				cluvec[indexDCRcol2]->colParticle(),
				cluvec[indexDCRacol1]->antiColParticle(),
				cluvec[indexDCRacol2]->antiColParticle()
				));

	cluvec[indexDCRcol1]->colParticle()->addChild(newClusterDCR);
	cluvec[indexDCRcol2]->colParticle()->addChild(newClusterDCR);
	cluvec[indexDCRacol1]->antiColParticle()->addChild(newClusterDCR);
	cluvec[indexDCRacol2]->antiColParticle()->addChild(newClusterDCR);

	// set new vertex
	newClusterDCR->setVertex(0.25*(
				cluvec[indexDCRcol1]->colParticle()->vertex() +
				cluvec[indexDCRcol2]->colParticle()->vertex() +
				cluvec[indexDCRacol1]->antiColParticle()->vertex() +
				cluvec[indexDCRacol2]->antiColParticle()->vertex()
				));

	// set beam remnants for new cluster
	if (cluvec[indexDCRcol1]->isBeamRemnant(colIdx[indexDCRcol1])) newClusterDCR->setBeamRemnant(0, true);
	if (cluvec[indexDCRcol2]->isBeamRemnant(colIdx[indexDCRcol2])) newClusterDCR->setBeamRemnant(1, true);
	if (cluvec[indexDCRacol1]->isBeamRemnant(anticolIdx[indexDCRacol1])) newClusterDCR->setBeamRemnant(2, true);
	if (cluvec[indexDCRacol2]->isBeamRemnant(anticolIdx[indexDCRacol2])) newClusterDCR->setBeamRemnant(3, true);

	return std::tuple <ClusterPtr, ClusterPtr> (newClusterMCR,newClusterDCR);
}

std::tuple  <ClusterPtr, ClusterPtr, ClusterPtr> ColourReconnector::_reconnect3Mto3M(
    ClusterPtr & c1, ClusterPtr & c2, ClusterPtr & c3, const int infos [3]) const {
  // check if mesonic clusters
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==2);
  assert(c3->numComponents()==2);

  ClusterVector oldclusters = {c1, c2, c3};
  ClusterVector newclusters;

  for (int i=0; i<3; i++) {

		if ( i!=infos[i] && _isColour8(oldclusters[i]->colParticle(),oldclusters[infos[i]]->antiColParticle())){
			std::cout << "Should never execute! _isColour8 in _reconnect3Mto3M "<< i <<" "<< infos[i]<<" check this earlier" << std::endl;
		}
		int c1_col=-1;
    int c2_anticol=-1;

    // get which index is coloured and which anticolour
    for (unsigned int ix=0; ix<2; ++ix) {
      if (oldclusters[i]->particle(ix)->hasColour(false)) c1_col  = ix;
      if (oldclusters[infos[i]]->particle(ix)->hasColour(true)) c2_anticol  = ix;
    }

    assert(c1_col>=0);
    assert(c2_anticol>=0);

    oldclusters[i]->colParticle()->abandonChild(oldclusters[i]);
    oldclusters[infos[i]]->antiColParticle()->abandonChild(oldclusters[infos[i]]);

    // form new cluster
    ClusterPtr newCluster = new_ptr(Cluster(oldclusters[i]->colParticle(), oldclusters[infos[i]]->antiColParticle()));

    oldclusters[i]->colParticle()->addChild(newCluster);
    oldclusters[infos[i]]->antiColParticle()->addChild(newCluster);

    // set new vertex
    newCluster->setVertex(0.5*(oldclusters[i]->colParticle()->vertex() +
                               oldclusters[infos[i]]->antiColParticle()->vertex()));

    // set beam remnants for new cluster
    if (oldclusters[i]->isBeamRemnant(c1_col)) newCluster->setBeamRemnant(0, true);
    if (oldclusters[infos[i]]->isBeamRemnant(c2_anticol)) newCluster->setBeamRemnant(1, true);
    newclusters.push_back(newCluster);
  }
  return std::tuple <ClusterPtr, ClusterPtr, ClusterPtr> (newclusters[0], newclusters[1], newclusters[2]);
}


pair  <ClusterPtr, ClusterPtr>
ColourReconnector::_reconnectMBtoMB(ClusterPtr & c1, ClusterPtr & c2, const int s0) const {
  // make c1 the mesonic cluster
  if (c1->numComponents()==2) {
    assert(c2->numComponents()==3);
  } else {
    return _reconnectMBtoMB(c2,c1,s0);
  }

  int c1_col=-1;
  int c1_anti=-1;
  // get which index is coloured and which anticolour
  for (unsigned int ix=0; ix<2; ++ix) {
    if (c1->particle(ix)->hasColour(false)) c1_col  = ix;
    else if (c1->particle(ix)->hasColour(true)) c1_anti = ix;

  }
  assert(c1_col>=0);
  assert(c1_anti>=0);

  // pointers for the new clusters
  ClusterPtr newCluster1;
  ClusterPtr newCluster2;
  if (c2->particle(0)->hasColour()==true) {
    // first case: we have a baryonic clusters

    // first make the new mesonic cluster
    c1->antiColParticle()->abandonChild(c1);
    c2->particleB(s0)->abandonChild(c2);

    newCluster1 = new_ptr(Cluster(c1->antiColParticle(), c2->particleB(s0)));

    c1->antiColParticle()->addChild(newCluster1);
    c2->particleB(s0)->addChild(newCluster1);

    // set new vertex
    newCluster1->setVertex(0.5*(c1->antiColParticle()->vertex() +
                                c2->particleB(s0)->vertex()));

    // set beam remnants for new cluster
    if (c1->isBeamRemnant(c1_anti)) newCluster1->setBeamRemnant(0, true);
    if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

    // then the baryonic one
    c1->colParticle()->abandonChild(c1);
    c2->particleB((s0+1)%3)->abandonChild(c2);
    c2->particleB((s0+2)%3)->abandonChild(c2);

    newCluster2 = new_ptr(Cluster(c1->colParticle(), c2->particleB((s0+1)%3), c2->particleB((s0+2)%3)));

    c1->colParticle()->addChild(newCluster2);
    c2->particleB((s0+1)%3)->addChild(newCluster2);
    c2->particleB((s0+2)%3)->addChild(newCluster2);

    // set new vertex
    newCluster2->setVertex(LorentzPoint());
  } else {
    // second case we have an antibaryonic cluster

    // first make the new mesonic cluster
    c1->colParticle()->abandonChild(c1);
    c2->particleB(s0)->abandonChild(c2);

    newCluster1 = new_ptr(Cluster(c1->colParticle(), c2->particleB(s0)));

    c1->colParticle()->addChild(newCluster1);
    c2->particleB(s0)->addChild(newCluster1);

    // set new vertex
    newCluster1->setVertex(0.5*(c1->colParticle()->vertex() +
                                c2->particleB(s0)->vertex()));

    // set beam remnants for new cluster
    if (c1->isBeamRemnant(c1_col)) newCluster1->setBeamRemnant(0, true);
    if (c2->isBeamRemnant(s0)) newCluster1->setBeamRemnant(1, true);

    // then the baryonic one
    c1->antiColParticle()->abandonChild(c1);
    c2->particleB((s0+1)%3)->abandonChild(c2);
    c2->particleB((s0+2)%3)->abandonChild(c2);

    newCluster2 =  new_ptr(Cluster(c1->antiColParticle(), c2->particleB((s0+1)%3), c2->particleB((s0+2)%3)));

    c1->antiColParticle()->addChild(newCluster2);
    c2->particleB((s0+1)%3)->addChild(newCluster2);
    c2->particleB((s0+2)%3)->addChild(newCluster2);

    // set new vertex
    newCluster2->setVertex(LorentzPoint());
  }
  return pair <ClusterPtr, ClusterPtr> (newCluster1, newCluster2);
}

void ColourReconnector::_2Bto2BreconnectionFinder(ClusterPtr & c1, ClusterPtr & c2,
    int & bswap1, int & bswap2, double min_displ_sum, string & kind_of_reco) const {
  double tmp_delta;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      // try swapping particle i of c1 with particle j of c2
      tmp_delta  = _displacementBaryonic(c2->particle(j), c1->particle((i+1)%3), c1->particle((i+2)%3));
      tmp_delta += _displacementBaryonic(c1->particle(i), c2->particle((j+1)%3), c2->particle((j+2)%3));

      if (tmp_delta < min_displ_sum) {
        // if minimal displacement select the 2Bto2B CR option
        min_displ_sum = tmp_delta;
        bswap1 = i;
        bswap2 = j;
        kind_of_reco = "2Bto2B";
      }
    }
  }

}

void ColourReconnector::_BbarBto3MreconnectionFinder(ClusterPtr & c1, ClusterPtr & c2,
    int & mswap0, int & mswap1, int & mswap2,
    double min_displ_sum, string & kind_of_reco) const {
  double pre_tmp_delta;
  double tmp_delta;
  for (int p1=0; p1 <3; p1++) {
    // make sure not to form a mesonic octet
    if (_isColour8(c1->particle(0), c2->particle(p1))) continue;

    pre_tmp_delta = _displacement(c1->particle(0), c2->particle(p1));
    for (int p2=1; p2<3; p2++) {

      // make sure not to form a mesonic octet
      if (_isColour8(c1->particle(1), c2->particle((p1+p2)%3))) continue;
      if (_isColour8(c1->particle(2), c2->particle(3-p1-((p1+p2)%3)))) continue;

      tmp_delta  = pre_tmp_delta + _displacement(c1->particle(1), c2->particle((p1+p2)%3));
      tmp_delta += 				 _displacement(c1->particle(2), c2->particle(3-p1-((p1+p2)%3)));

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_delta *=_mesonToBaryonFactor;

      if (tmp_delta < min_displ_sum) {
        // if minimal displacement select the 2Bto3M CR option
        min_displ_sum = tmp_delta;
        mswap0 = p1;
        mswap1 = (p1+p2)%3;
        mswap2 = 3-p1-((p1+p2)%3);
        kind_of_reco = "2Bto3M";

      }
    }
  }
}

void ColourReconnector::_BMtoBMreconnectionfinder(ClusterPtr & c1, ClusterPtr & c2, int & swap, double min_displ_sum,
    string & kind_of_reco) const {
  assert(c1->numComponents()==2);
  assert(c2->numComponents()==3);
  double tmp_displ = 0;
  for (int i=0; i<3; i++) {
    // Differ if the second cluster is baryonic or antibaryonic
    if (c2->particle(0)->hasColour()) {
      // c2 is baryonic

      // veto mesonic octets
      if (_isColour8(c2->particle(i), c1->antiColParticle())) continue;

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_displ = _mesonToBaryonFactor * _displacement(c2->particle(i), c1->antiColParticle());
      tmp_displ += _displacementBaryonic(c1->colParticle(), c2->particle((i+1)%3), c2->particle((i+2)%3));
    } else {
      // c2 is antibaryonic

      // veto mesonic octets
      if (_isColour8(c2->particle(i), c1->colParticle())) continue;

      // factor _mesonToBaryonFactor to compare Baryonic an mesonic cluster
      tmp_displ = _mesonToBaryonFactor * _displacement(c2->particle(i), c1->colParticle());
      tmp_displ *= _displacementBaryonic(c1->antiColParticle(), c2->particle((i+1)%3), c2->particle((i+2)%3));
    }
    if (tmp_displ < min_displ_sum) {
      // if minimal displacement select the MBtoMB CR option
      min_displ_sum = tmp_displ;
      swap = i;
      kind_of_reco = "MBtoMB";
    }
  }
  return;
}

void ColourReconnector::_3MtoXreconnectionfinder(std::vector<CluVecIt> & cv, int & swap0, int & swap1,
    int & swap2, double min_displ_sum, string & kind_of_reco) const {
  // case of 3M->BbarB CR
  double _tmp_displ;
  _tmp_displ  = _displacementBaryonic((*cv[0])->colParticle(),     (*cv[1])->colParticle(),     (*cv[2])->colParticle());
  _tmp_displ += _displacementBaryonic((*cv[0])->antiColParticle(), (*cv[1])->antiColParticle(), (*cv[2])->antiColParticle());
  if (_tmp_displ < min_displ_sum) {
    // if minimal displacement select the 3Mto2B CR option
    kind_of_reco = "3Mto2B";
    min_displ_sum = _tmp_displ;
  }
  // case for 3M->3M CR
  /**
   * if 3Mto3M reco probability (_preco3M_3M) is 0 we skip this loop
   * since no 3Mto3M CR shall be performed
   */
  int i,j;
  int i1,i2,i3;
  for (i = 0; _preco3M_3M && i<3; i++) {
    // veto mesonic octets
    if (_isColour8((*cv[0])->colParticle(), (*cv[i])->antiColParticle())) continue;

    // factor _mesonToBaryonFactor to compare baryonic an mesonic cluster
    _tmp_displ = _mesonToBaryonFactor * _displacement((*cv[0])->colParticle(), (*cv[i])->antiColParticle());
    for (j=1; j<3; j++) {
      // i1, i2, i3 are pairwise distinct
      i1=i;
      i2=((j+i)%3);
      if (i1==0 && i2==1) continue;
      i3=(3-i-((j+i)%3));

      // veto mesonic octets
      if (_isColour8((*cv[1])->colParticle(), (*cv[i2])->antiColParticle())) continue;
      if (_isColour8((*cv[2])->colParticle(), (*cv[i3])->antiColParticle())) continue;

      _tmp_displ += _mesonToBaryonFactor * _displacement((*cv[1])->colParticle(), (*cv[i2])->antiColParticle());
      _tmp_displ += _mesonToBaryonFactor * _displacement((*cv[2])->colParticle(), (*cv[i3])->antiColParticle());

      if (_tmp_displ < min_displ_sum) {
        // if minimal displacement select the 3Mto3M CR option
        kind_of_reco = "3Mto3M";
        min_displ_sum = _tmp_displ;
        swap0 = i1;
        swap1 = i2;
        swap2 = i3;
      }
    }
  }
}


pair <int,int> ColourReconnector::_shuffle(
    const PVector & q, const PVector & aq, unsigned maxtries) const {

  const size_t nclusters = q.size();
  assert (nclusters > 1);
  assert (aq.size() == nclusters);

  int i, j;
  unsigned tries = 0;
  bool octet=false;

  do {
    // find two different random integers in the range [0, nclusters)
    i = UseRandom::irnd( nclusters );
    do {
      j = UseRandom::irnd( nclusters );
    } while (i == j);

    // check if one of the two potential clusters would be a colour octet state
    octet = _isColour8( q[i], aq[j] ) || _isColour8( q[j], aq[i] ) ;
    tries++;
  } while (octet && tries < maxtries);

  if (octet) i = j = -1;
  return make_pair(i,j);
}



bool ColourReconnector::_becomesColour8Cluster(const ClusterPtr & c1, const ClusterPtr & c2) const {
	return _isColour8(c1->colParticle(),c2->antiColParticle()) || _isColour8(c2->colParticle(),c1->antiColParticle());
}
bool ColourReconnector::_isColour8(tcPPtr p, tcPPtr q) const {
  bool octet = false;
	if(_octetOption<0) return octet;

  // make sure we have a triplet and an anti-triplet
  if ( ( p->hasColour() && q->hasAntiColour() ) ||
       ( p->hasAntiColour() && q->hasColour() ) ) {

    // true if p and q are originated from a colour octet
    if ( !p->parents().empty() && !q->parents().empty() ) {
      octet = ( p->parents()[0] == q->parents()[0] ) &&
              ( p->parents()[0]->data().iColour() == PDT::Colour8 );
    }

    // (Final) option: check if same colour8 parent
    // or already found an octet.
    if(_octetOption==0||octet) return octet;

    // (All) option handling more octets
    // by browsing particle history/colour lines.
    tColinePtr cline,aline;

    // Get colourlines form final states.
    if(p->hasColour() && q->hasAntiColour()) {
      cline  = p->    colourLine();
      aline  = q->antiColourLine();
    } else {
      cline  = q->    colourLine();
      aline  = p->antiColourLine();
    }

    // Follow the colourline of p.
    if ( !p->parents().empty() ) {
      tPPtr parent = p->parents()[0];
      while (parent) {
        if(parent->data().iColour() == PDT::Colour8) {
          // Coulour8 particles should have a colour
          // and an anticolour line. Currently the
          // remnant has none of those. Since the children
          // of the remnant are not allowed to emit currently,
          // the colour octet remnant is handled by the return
          // statement above. The assert also catches other
          // colour octets without clines. If the children of
          // a remnant should be allowed to emit, the remnant
          // should get appropriate colour lines and
          // colour states.
          //  See Ticket: #407
          //  assert(parent->colourLine()&&parent->antiColourLine());
          octet = (parent->    colourLine()==cline &&
                   parent->antiColourLine()==aline);
        }
        if(octet||parent->parents().empty()) break;
        parent = parent->parents()[0];
      }
    }
  }

  return octet;
}


void ColourReconnector::persistentOutput(PersistentOStream & os) const {
  os
      << _hadronSpectrum
      << _clusterFinder
      << _clreco
      << _crIterations
      << _algorithm
      << _annealingFactor
      << _annealingSteps
      << _triesPerStepFactor
      << _initTemp
      << _precoMesonic
      << _precoBaryonic
      << _precoDiquark
      << _preco3M_3M
      << _preco3M_BBbar
      << _precoBbarB_3M
      << _preco2B_2B
      << _precoMB_MB
      << _stepFactor
      << _mesonToBaryonFactor
      << ounit(_maxDistance, femtometer)
      << _octetOption
      << _localCR
      << _causalCR
      << _debug
      << _junctionMBCR
			<< _dynamicCR
			<< _diquarkCR
			<< _dynamicCRscale
			<< _dynamicCRalphaS
			<< _phaseSpaceDiquarkFission
			<< _cutDiquarkClusterFormation
      ;
}

void ColourReconnector::persistentInput(PersistentIStream & is, int) {
  is
      >> _hadronSpectrum
      >> _clusterFinder
      >> _clreco
      >> _crIterations
      >> _algorithm
      >> _annealingFactor
      >> _annealingSteps
      >> _triesPerStepFactor
      >> _initTemp
      >> _precoMesonic
      >> _precoBaryonic
      >> _precoDiquark
      >> _preco3M_3M
      >> _preco3M_BBbar
      >> _precoBbarB_3M
      >> _preco2B_2B
      >> _precoMB_MB
      >> _stepFactor
      >> _mesonToBaryonFactor
      >> iunit(_maxDistance, femtometer)
      >> _octetOption
      >> _localCR
      >> _causalCR
      >> _debug
      >> _junctionMBCR
			>> _dynamicCR
			>> _diquarkCR
			>> _dynamicCRscale
			>> _dynamicCRalphaS
			>> _phaseSpaceDiquarkFission
			>> _cutDiquarkClusterFormation
      ;
}


void ColourReconnector::Init() {
  static ClassDocumentation<ColourReconnector> documentation
  ("This class is responsible of the colour reconnection.");

  // Reference to the HadronSpectrum
  static Reference<ColourReconnector,HadronSpectrum>
    interfaceHadronSpectrum("HadronSpectrum",
		      "A reference to the object HadronSpectrum"
					" Needed for thresholds",
		      &Herwig::ColourReconnector::_hadronSpectrum,
		      false, false, true, false);

  // Reference to the ClusterFinder
  static Reference<ColourReconnector,ClusterFinder>
    interfaceClusterFinder("ClusterFinder",
		      "A reference to the object ClusterFinder "
					"Needed for Diquark On shell treatment",
		      &Herwig::ColourReconnector::_clusterFinder,
		      false, false, true, false);

  // Global switch for switching off CR alltogheter
  static Switch<ColourReconnector,int> interfaceColourReconnection
  ("ColourReconnection",
   "Colour reconnections",
   &ColourReconnector::_clreco, 0, true, false);
  static SwitchOption interfaceColourReconnectionNo
  (interfaceColourReconnection,
   "No",
   "Colour reconnections off",
   0);
  static SwitchOption interfaceColourReconnectionYes
  (interfaceColourReconnection,
   "Yes",
   "Colour reconnections on",
   1);

  // Number of iterations for multiple pass-throughs of the CR algorithm
  static Parameter<ColourReconnector, unsigned int> interfaceColourReconnectionIterations
  ("ColourReconnectionIterations",
   "Choose the number of iterations the chosen CR algorithm is performed",
   &ColourReconnector::_crIterations, 1, 1, 1, 100,
   false, false, Interface::limited);

  // Algorithm interface
  static Switch<ColourReconnector, int> interfaceAlgorithm
  ("Algorithm",
   "Specifies the colour reconnection algorithm",
   &ColourReconnector::_algorithm, 0, true, false);
  static SwitchOption interfaceAlgorithmPlain
  (interfaceAlgorithm,
   "Plain",
   "Plain colour reconnection as in Herwig 2.5.0",
   0);
  static SwitchOption interfaceAlgorithmStatistical
  (interfaceAlgorithm,
   "Statistical",
   "Statistical colour reconnection using simulated annealing",
   1);
  static SwitchOption interfaceAlgorithmBaryonic
  (interfaceAlgorithm,
   "Baryonic",
   "Baryonic cluster reconnection",
   2);
  static SwitchOption interfaceAlgorithmBaryonicMesonic
  (interfaceAlgorithm,
   "BaryonicMesonic",
   "Baryonic cluster reconnection with reconnections to and from Mesonic Clusters",
   3);
  static SwitchOption interfaceAlgorithmBaryonicDiquarkCluster
  (interfaceAlgorithm,
   "BaryonicDiquarkCluster",
   "Baryonic colour reconnection which allows for the formation of DiquarkCluster-like CR",
   4);

  // Turn on/off the generation of diquark cluster via CR
  static Switch<ColourReconnector,int> interfaceColourDiquarkCR
  ("DiquarkCR",
   "Allow diquark type colour Reconnection."
	 "NOTE: Necessary to be Yes if BaryonicDiquarkCluster algorithm is chosen.",
   &ColourReconnector::_diquarkCR, 0, true, false);
  static SwitchOption interfaceDiquarkCRNo
  (interfaceColourDiquarkCR,
   "No",
   "Forbid diquark type colour reconnections",
   0);
  static SwitchOption interfaceDiquarkCRYes
  (interfaceColourDiquarkCR,
   "Yes",
   "Allow diquark type colour reconnections",
   1);

  static Switch<ColourReconnector,int> interfaceColourDynamicCR
  ("DynamicCR",
   "Use dynamic weight for Colour reconnections defined by soft gluon evolution"
	 "\nNOTE: Only availible for Plain, Baryonic, BaryonicDiquarkCluster algorithm.",
   &ColourReconnector::_dynamicCR, 0, true, false);
  static SwitchOption interfaceDynamicCRNo
  (interfaceColourDynamicCR,
   "No",
   "Use regular CR with fixed probabilities",
   0);
  static SwitchOption interfaceDynamicCRYes
  (interfaceColourDynamicCR,
   "Yes",
   "Use dynamic CR with kinematic dependent probabilities",
   1);
  static SwitchOption interfaceDynamicCRNotExponentiated
  (interfaceColourDynamicCR,
   "NotExponentiated",
   "Use dynamic CR with kinematic dependent probabilities"
	 ", but without exponentiated soft anomalous dimension."
	 "NOTE: Only for testing.",
   2);

  // General Parameters and switches
  static Parameter<ColourReconnector, double> interfaceDynamicScale
  ("DynamicScale",
   "Choose dynamic scale of soft gluon evolution for DynamicCR where"
	 " mu = DynamicScale*(mConstU+mConstU)",
   &ColourReconnector::_dynamicCRscale, 1.0, 1e-14, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceDynamicAlphaS
  ("DynamicAlphaS",
   "Choose dynamic alphaS of soft gluon evolution for DynamicCR",
   &ColourReconnector::_dynamicCRalphaS, 0.8, 0.001, 10.0,
   false, false, Interface::limited);


  // Statistical CR Parameters:
  static Parameter<ColourReconnector, double> interfaceMtrpAnnealingFactor
  ("AnnealingFactor",
   "The annealing factor is the ratio of the temperatures in two successive "
   "temperature steps.",
   &ColourReconnector::_annealingFactor, 0.9, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,unsigned> interfaceMtrpAnnealingSteps
  ("AnnealingSteps",
   "Number of temperature steps in the statistical annealing algorithm",
   &ColourReconnector::_annealingSteps, 50, 1, 10000,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,double> interfaceMtrpTriesPerStepFactor
  ("TriesPerStepFactor",
   "The number of reconnection tries per temperature steps is the number of "
   "clusters times this factor.",
   &ColourReconnector::_triesPerStepFactor, 5.0, 0.0, 100.0,
   false, false, Interface::limited);


  static Parameter<ColourReconnector,double> interfaceMtrpInitialTemp
  ("InitialTemperature",
   "Factor used to determine the initial temperature from the median of the "
   "energy change in a few random rearrangements.",
   &ColourReconnector::_initTemp, 0.1, 0.00001, 100.0,
   false, false, Interface::limited);




  // Plain and Baryonic CR Paramters
  static Parameter<ColourReconnector, double> interfaceRecoProb
  ("ReconnectionProbability",
   "Probability that a found two meson to two meson reconnection possibility is actually accepted (used in Plain & Baryonic)",
   &ColourReconnector::_precoMesonic, 0.5, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector,double> interfaceRecoProbBaryonic
  ("ReconnectionProbabilityBaryonic",
   "Probability that a found reconnection possibility is actually accepted (used in Baryonic)",
   &ColourReconnector::_precoBaryonic, 0.5, 0.0, 1.0,
   false, false, Interface::limited);


  static Parameter<ColourReconnector,double> interfaceRecoProbDiquark
  ("ReconnectionProbabilityDiquark",
   "Probability for forming a tetra-quark cluster",
   &ColourReconnector::_precoDiquark, 0.5, 0.0, 1.0,
   false, false, Interface::limited);


  // BaryonicMesonic CR Paramters
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability3Mto3M
  ("ReconnectionProbability3Mto3M",
   "Probability that a reconnection candidate is accepted for reconnecting 3M -> 3M\'",
   &ColourReconnector::_preco3M_3M, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability3MtoBBbar
  ("ReconnectionProbability3MtoBBbar",
   "Probability that a reconnection candidate is accepted for reconnecting 3M -> B,Bbar",
   &ColourReconnector::_preco3M_BBbar, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbabilityBbarBto3M
  ("ReconnectionProbabilityBbarBto3M",
   "Probability that a reconnection candidate is accepted for reconnecting B,Bbar -> 3M",
   &ColourReconnector::_precoBbarB_3M, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbability2Bto2B
  ("ReconnectionProbability2Bto2B",
   "Probability that a reconnection candidate is accepted for reconnecting 2B -> 2B\' or 2Bbar -> 2Bbar\'",
   &ColourReconnector::_preco2B_2B, 0.5, 0.0, 1.0,
   false, false, Interface::limited);
  static Parameter<ColourReconnector, double> interfaceReconnectionProbabilityMBtoMB
  ("ReconnectionProbabilityMBtoMB",
   "Probability that a reconnection candidate is accepted for reconnecting M,B -> M\',B\' or M,Bbar -> M\',Bbar\'",
   &ColourReconnector::_precoMB_MB, 0.5, 0.0, 1.0,
   false, false, Interface::limited);

  static Parameter<ColourReconnector, double> interfaceFactorforStep
  ("StepFactor",
   "Factor for how many reconnection-tries are made in the BaryonicMesonic algorithm",
   &ColourReconnector::_stepFactor, 1.0, 0.11111, 10.,
   false, false, Interface::limited);// at least 3 Clusters -> _stepFactorMin=1/9

  static Parameter<ColourReconnector, double> interfaceMesonToBaryonFactor
  ("MesonToBaryonFactor",
   "Factor for comparing mesonic clusters to baryonic clusters in the displacement if BaryonicMesonic CR model is chosen",
   &ColourReconnector::_mesonToBaryonFactor, 2.0, 0.01, 100.0,
   false, false, Interface::limited);


  // General Parameters and switches
  static Parameter<ColourReconnector, Length> interfaceMaxDistance
  ("MaxDistance",
   "Maximum distance between the clusters at which to consider rearrangement"
   " to avoid colour reconneections of displaced vertices (used in all Algorithms). No unit means femtometer",
   &ColourReconnector::_maxDistance, femtometer, 1000.*femtometer, 0.0*femtometer, 1e100*femtometer,
   false, false, Interface::limited);


  static Switch<ColourReconnector, int> interfaceOctetTreatment
  ("OctetTreatment",
   "Which octets are not allowed to be reconnected (used in all Algorithms)",
   &ColourReconnector::_octetOption, 0, false, false);
  static SwitchOption interfaceOctetTreatmentFinal
  (interfaceOctetTreatment,
   "Final",
   "Only prevent for the final (usuaslly non-perturbative) g -> q qbar splitting",
   0);
  static SwitchOption interfaceOctetTreatmentAll
  (interfaceOctetTreatment,
   "All",
   "Prevent for all octets",
   1);
  static SwitchOption interfaceOctetTreatmentNone
  (interfaceOctetTreatment,
   "None",
   "Accept all octets. "
	 "NOTE: If a static gluon constituent mass is chosen this option is unphysical. "
	 "It will lead to mConstGluon fixed mass clusters!",
   -1);

  static Switch<ColourReconnector, int> interfaceLocalCR
  ("LocalCR",
   "Option for colour reconnecting only if clusters are less distant than MaxDistance",
   &ColourReconnector::_localCR, 0, true, false);
  static SwitchOption interfaceLocalCRYes
  (interfaceLocalCR,
   "Yes",
   "activate spatial veto",
   1);
  static SwitchOption interfaceLocalCRNo
  (interfaceLocalCR,
   "No",
   "deactivate spatial veto",
   0);
  static Switch<ColourReconnector, int> interfaceCausalCR
  ("CausalCR",
   "Option for colour reconnecting only if clusters their vertices "
   "have a positive spacetime difference",
   &ColourReconnector::_causalCR, 0, true, false);
  static SwitchOption interfaceCausalCRYes
  (interfaceCausalCR,
   "Yes",
   "enable causal veto",
   1);
  static SwitchOption interfaceCausalCRNo
  (interfaceCausalCR,
   "No",
   "disable causal veto",
   0);

  static Switch<ColourReconnector, int> interfaceJunction
  ("Junction",
   "Option for using Junction-like displacement in rapidity-phi "
   "plane to compare baryonic cluster "
   "instead of pairwise distance (for BaryonicMesonic model)",
   &ColourReconnector::_junctionMBCR, 1, true, false);
  static SwitchOption interfaceJunctionYes
  (interfaceJunction,
   "Yes",
   "Using junction-like model instead of pairwise distance model",
   1);
  static SwitchOption interfaceJunctionNo
  (interfaceJunction,
   "No",
   "Using pairwise distance model instead of junction-like model",
   0);

  // Debug
  static Switch<ColourReconnector, int> interfaceDebug
  ("Debug",
   "Make a file with some Information of the BaryonicMesonic Algorithm",
   &ColourReconnector::_debug, 0, true, false);
  static SwitchOption interfaceDebugNo
  (interfaceDebug,
   "No",
   "Debug Information for ColourReconnector Off",
   0);
  static SwitchOption interfaceDebugYes
  (interfaceDebug,
   "Yes",
   "Debug Information for ColourReconnector On",
   1);

  static Switch<ColourReconnector,int> interfacePhaseSpaceDiquarkFission
  ("PhaseSpaceDiquarkFission",
   "Only for dynamic colour reconnection choose if capturing cluster decay"
	 " phase space for formation of a diquark cluster in the transition probabilities",
   &ColourReconnector::_phaseSpaceDiquarkFission, 0, true, false);
  static SwitchOption interfacePhaseSpaceDiquarkFissionNo
  (interfacePhaseSpaceDiquarkFission,
   "No",
   "Not adding the decay phasespace to Diquark Colour Reconnection",
   0);
  static SwitchOption interfacePhaseSpaceDiquarkFissionYes
  (interfacePhaseSpaceDiquarkFission,
   "Yes",
   "Adding the decay phasespace to Diquark Colour Reconnection",
   1);
  static SwitchOption interfacePhaseSpaceDiquarkFissionConstituentMasses
  (interfacePhaseSpaceDiquarkFission,
   "ConstituentMasses",
   "Adding the decay phasespace to Diquark Colour Reconnection",
   2);

  static Parameter<ColourReconnector, double> interfaceCutDiquarkClusterFormation
  ("CutDiquarkClusterFormation",
   "Cut on diquark cluster formation such that clusters are accepted only if "
	 "(p1*p2/(m1*m2)-1)<cut for diquark and antidiquark. Note that setting this to 0"
	 " applies no cut!",
   &ColourReconnector::_cutDiquarkClusterFormation, 0.5, 0.0, 100.0,
   false, false, Interface::limited);

  static Switch<ColourReconnector,bool> interfaceDoublyHeavyBaryon
    ("DoublyHeavyBaryon",
     "Allow the production of doubly heavy baryons in colour reconnection",
     &ColourReconnector::_doublyHeavyBaryon, true, false, false);
  static SwitchOption interfaceDoublyHeavyBaryonYes
    (interfaceDoublyHeavyBaryon,
     "Yes",
     "Allow doubly heavy baryon production",
     true);
  static SwitchOption interfaceDoublyHeavyBaryonNo
    (interfaceDoublyHeavyBaryon,
     "No",
     "Don't allow doubly heavy baryon production",
     false);

}

