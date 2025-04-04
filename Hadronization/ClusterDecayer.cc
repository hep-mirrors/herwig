// -*- C++ -*-
//
// ClusterDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterDecayer class.
//

#include "ClusterDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Herwig/Utilities/Kinematics.h"
#include "Cluster.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include <ThePEG/Repository/UseRandom.h>
#include "ThePEG/Interface/ParMap.h"

#include "Herwig/Utilities/Histogram.h"

using namespace Herwig;

DescribeClass<ClusterDecayer,Interfaced>
describeClusterDecayer("Herwig::ClusterDecayer","Herwig.so");

ClusterDecayer::ClusterDecayer() :
  _clDirLight(1),
  _clDirExotic(1),
  _clSmrLight(0.0),
  _clSmrExotic(0.0),
  _onshell(false),
  _masstry(20),
  _kinematics(0)
{}

IBPtr ClusterDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterDecayer::fullclone() const {
  return new_ptr(*this);
}

void ClusterDecayer::persistentOutput(PersistentOStream & os) const
{
  os << _clDirLight << _clDirHeavy
     << _clDirExotic << _clSmrLight << _clSmrHeavy
     << _clSmrExotic << _onshell << _masstry
		 << _kinematics
     << _hadronSpectrum;
}

void ClusterDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _clDirLight >> _clDirHeavy
     >> _clDirExotic >> _clSmrLight >> _clSmrHeavy
     >> _clSmrExotic >> _onshell >> _masstry
		 >> _kinematics
     >> _hadronSpectrum;
}

void ClusterDecayer::doinit() {
  Interfaced::doinit();
  for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
    if ( _clSmrHeavy.find(id) == _clSmrHeavy.end() ||
	 _clDirHeavy.find(id) == _clDirHeavy.end() )
      throw InitException() << "not all parameters have been set for heavy quark cluster decay";
  }
}

void ClusterDecayer::Init() {

  static ClassDocumentation<ClusterDecayer> documentation
    ("This class is responsible for the two-body decays of normal clusters");

  static Reference<ClusterDecayer,HadronSpectrum> interfaceHadronSpectrum
    ("HadronSpectrum",
     "Set the hadron spectrum for this parton splitter.",
     &ClusterDecayer::_hadronSpectrum, false, false, true, false);

  //ClDir for light, Bottom, Charm and exotic (e.g Susy) quarks

  static Switch<ClusterDecayer,bool> interfaceClDirLight
    ("ClDirLight",
     "Cluster direction for light quarks",
     &ClusterDecayer::_clDirLight, true, false, false);
  static SwitchOption interfaceClDirLightPreserve
    (interfaceClDirLight,
     "Preserve",
     "Preserve the direction of the quark from the perturbative stage"
     " as the direction of the meson produced from it",
     true);
  static SwitchOption interfaceClDirLightIsotropic
    (interfaceClDirLight,
     "Isotropic",
     "Assign the direction of the meson randomly",
     false);
  
  static ParMap<ClusterDecayer,bool> interfaceClDirHeavy
    ("ClDirHeavy",
     "Preserve the direction of the given heavy quark from the perturbative stage.",
     &ClusterDecayer::_clDirHeavy, -1, true, false, true,
     false, false, Interface::limited);

  static Switch<ClusterDecayer,bool> interfaceClDirExotic
    ("ClDirExotic",
     "Cluster direction for exotic quarks",
     &ClusterDecayer::_clDirExotic, true, false, false);
  static SwitchOption interfaceClDirExoticPreserve
    (interfaceClDirExotic,
     "Preserve",
     "Preserve the direction of the quark from the perturbative stage"
     " as the direction of the meson produced from it",
     true);
  static SwitchOption interfaceClDirExoticIsotropic
    (interfaceClDirExotic,
     "Isotropic",
     "Assign the direction of the meson randomly",
     false);


  static Switch<ClusterDecayer,int> interfaceKinematics
    ("Kinematics",
     "Choose kinematics of Cluster Decay",
     &ClusterDecayer::_kinematics, 0, true, false);
  static SwitchOption interfaceKinematicsDefault
    (interfaceKinematics,
     "Default",
     "default kinematics. Isotropic for soft Clusters"
		 ", smeared aligned for perturbative Clusters",
     0);
  static SwitchOption interfaceKinematicsTchannel
    (interfaceKinematics,
     "Tchannel",
     "Experimental!: t-channel like kinematics "
		 "using the matrix element 1/((p1.p3)^2 - m1*m3)^2",
     1);
  static SwitchOption interfaceKinematicsAllIsotropic
    (interfaceKinematics,
     "AllIsotropic",
		 "All clusters decay isotropically"
     "NOTE: Just for Debug and Testing",
     2);
  static SwitchOption interfaceKinematicsAllAligned
    (interfaceKinematics,
     "AllAligned",
     "All clusters decay in the direction of their constituents"
     "NOTE: Just for Debug and Testing",
     3);
  static SwitchOption interfaceKinematicsTchannelHemisphere
    (interfaceKinematics,
     "TchannelHemisphere",
     "Experimental!: t-channel like kinematics "
		 "using the matrix element 1/((p1.p3)^2 - m1*m3)^2"
		 "But limited the sampling to the hemisphere of the original constituents"
     "NOTE: Just for Debug and Testing",
     -1);
  static SwitchOption interfaceKinematicsAllIsotropicHemisphere
    (interfaceKinematics,
     "IsotropicHemisphere",
     "All clusters decay semi-isotropically "
		 "limited the sampling to the hemisphere of the original constituents"
     "NOTE: Just for Debug and Testing",
     -2);

  // ClSmr for ligth, Bottom, Charm and exotic (e.g. Susy) quarks
  static Parameter<ClusterDecayer,double>
    interfaceClSmrLight ("ClSmrLight", "cluster direction Gaussian smearing for light quark",
                     &ClusterDecayer::_clSmrLight, 0, 0.0 , 0.0 , 2.0,false,false,false);
  static ParMap<ClusterDecayer,double> interfaceClSmrHeavy
    ("ClSmrHeavy",
     "cluster direction Gaussian smearing for heavy quarks",
     &ClusterDecayer::_clSmrHeavy, -1, 0.0, 0.0, 2.0,
     false, false, Interface::limited);
  static Parameter<ClusterDecayer,double>
    interfaceClSmrExotic ("ClSmrExotic", "cluster direction Gaussian smearing for exotic quark",
                     &ClusterDecayer::_clSmrExotic, 0, 0.0 , 0.0 , 2.0,false,false,false);

  static Switch<ClusterDecayer,bool> interfaceOnShell
    ("OnShell",
     "Whether or not the hadrons produced should by on shell or generated using the"
     " mass generator.",
     &ClusterDecayer::_onshell, false, false, false);
  static SwitchOption interfaceOnShellOnShell
    (interfaceOnShell,
     "Yes",
     "Produce the hadrons on shell",
     true);
  static SwitchOption interfaceOnShellOffShell
    (interfaceOnShell,
     "No",
     "Generate the masses using the mass generator.",
     false);

  static Parameter<ClusterDecayer,unsigned int> interfaceMassTry
    ("MassTry",
     "The number attempts to generate the masses of the hadrons produced"
     " in the cluster decay.",
     &ClusterDecayer::_masstry, 20, 1, 50,
     false, false, Interface::limited);

}


void ClusterDecayer::decay(const ClusterVector & clusters, tPVector & finalhadrons) {
  // Loop over all clusters, and if they are not too heavy (that is
  // intermediate clusters that have undergone to fission) or not
  // too light (that is final clusters that have been already decayed
  // into single hadron) then decay them into two hadrons.

  for (ClusterVector::const_iterator it = clusters.begin();
	 it != clusters.end(); ++it) {
    if ((*it)->isAvailable() && !(*it)->isStatusFinal()
	&& (*it)->isReadyToDecay()) {
      pair<PPtr,PPtr> prod = decayIntoTwoHadrons(*it);
      if(!prod.first)
	throw Exception() << "Can't perform decay of cluster ("
			  << (*it)->particle(0)->dataPtr()->PDGName() << " "
			  << (*it)->particle(1)->dataPtr()->PDGName() << ")\nThreshold = " 
			  << spectrum()->massLightestHadronPair((*it)->particle(0)->dataPtr(),
                                                                (*it)->particle(1)->dataPtr())/GeV
			  << "\nLightest hadron pair = ( "
			  << spectrum()->lightestHadronPair((*it)->particle(0)->dataPtr(),
					  						  (*it)->particle(1)->dataPtr()).first->PDGName() << ",\t"
			  << spectrum()->lightestHadronPair((*it)->particle(0)->dataPtr(),
					  						  (*it)->particle(1)->dataPtr()).second->PDGName() << " )\nCluster =\n"

			  << **it
			  << "\nMass = " << (*it)->mass()/GeV
			  << "in ClusterDecayer::decay()"
			  << Exception::eventerror;
      (*it)->addChild(prod.first);
      (*it)->addChild(prod.second);
      finalhadrons.push_back(prod.first);
      finalhadrons.push_back(prod.second);
    }
  }
}


pair<PPtr,PPtr> ClusterDecayer::decayIntoTwoHadrons(tClusterPtr ptr) {
	switch (abs(_kinematics))
	{
		case 0:
			// Default
			return decayIntoTwoHadronsDefault(ptr);
		case 1:
			// Tchannel- like
			return decayIntoTwoHadronsNew(ptr);
		case 2:
			// AllIsotropic
			return decayIntoTwoHadronsNew(ptr);
		case 3:
			// AllAligned
			return decayIntoTwoHadronsNew(ptr);
		default:
			assert(false);
	}
}

pair<PPtr,PPtr> ClusterDecayer::decayIntoTwoHadronsDefault(tClusterPtr ptr) {
	using Constants::pi;
	using Constants::twopi;
	// To decay the cluster into two hadrons one must distinguish between
	// constituent quarks (or diquarks) that originate from perturbative
	// processes (hard process or parton shower) from those that are generated
	// by the non-perturbative gluon splitting or from fission of heavy clusters.
	// In the latter case the two body decay is assumed to be *isotropic*.
	// In the former case instead, if proper flag are activated, the two body
	// decay is assumed to "remember" the direction of the constituents inside
	// the cluster, in the cluster frame. The actual smearing of the hadron
	// directions around the direction of the constituents, in the cluster
	// frame, can be different between non-b hadrons and b-hadrons, but is given
	// by the same functional form:
	//          cosThetaSmearing = 1 + smearFactor * log( rnd() )
	// (repeated until cosThetaSmearing > -1)
	// where the factor smearFactor is different between b- and non-b hadrons.
	//
	// We need to require (at least at the moment, maybe in the future we
	// could change it) that the cluster has exactly two components.
	// If this is not the case, then send a warning because it is not suppose
	// to happen, and then return.
	if ( ptr->numComponents() != 2 ) {
		generator()->logWarning( Exception("ClusterDecayer::decayIntoTwoHadrons "
					"***Still cluster with not exactly 2 components*** ",
					Exception::warning) );
		return pair<PPtr,PPtr>();
	}

	// Extract the id and particle pointer of the two components of the cluster.
	tPPtr ptr1 = ptr->particle(0);
	tPPtr ptr2 = ptr->particle(1);
	tcPDPtr ptr1data = ptr1->dataPtr();
	tcPDPtr ptr2data = ptr2->dataPtr();

	bool isHad1FlavSpecial    = false;
	bool cluDirHad1      = _clDirLight;
	double cluSmearHad1 = _clSmrLight;
	bool isHad2FlavSpecial    = false;
	bool cluDirHad2      = _clDirLight;
	double cluSmearHad2 = _clSmrLight;

	if (spectrum()->isExotic(ptr1data)) {
		isHad1FlavSpecial  = true;
		cluDirHad1   = _clDirExotic;
		cluSmearHad1 = _clSmrExotic;
	} else {
		for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
			if ( spectrum()->hasHeavy(id,ptr1data) ) {
				cluDirHad1   = _clDirHeavy[id];
				cluSmearHad1 = _clSmrHeavy[id];
			}
		}
	}


	if (spectrum()->isExotic(ptr2data)) {
		isHad2FlavSpecial  = true;
		cluDirHad2   = _clDirExotic;
		cluSmearHad2 = _clSmrExotic;
	} else {
		for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
			if ( spectrum()->hasHeavy(id,ptr2data) ) {
				cluDirHad2   = _clDirHeavy[id];
				cluSmearHad2 = _clSmrHeavy[id];
			}
		}
	}

	bool isOrigin1Perturbative = ptr->isPerturbative(0);
	bool isOrigin2Perturbative = ptr->isPerturbative(1);

	// We have to decide which, if any, of the two hadrons will have
	// the momentum, in the cluster parent frame, smeared around the
	// direction of its constituent (for Had1 is the one pointed by
	// ptr1, and for Had2 is the one pointed by ptr2).
	// This happens only if the flag _ClDirX is 1 and the constituent is
	// perturbative (that is not coming from nonperturbative gluon splitting
	// or cluster fission). In the case that both the hadrons satisfy this
	// two requirements (of course only one must be treated, because the other
	// one will have the momentum automatically fixed by the momentum
	// conservation) then more priority is given in the case of a b-hadron.
	// Finally, in the case that the two hadrons have same priority, then
	// we choose randomly, with equal probability, one of the two.

	int priorityHad1 = 0;
	if ( cluDirHad1  &&  isOrigin1Perturbative ) {
		priorityHad1 = isHad1FlavSpecial ? 2 : 1;
	}
	int priorityHad2 = 0;
	if ( cluDirHad2  &&  isOrigin2Perturbative ) {
		priorityHad2 = isHad2FlavSpecial ? 2 : 1;
	}
	if ( priorityHad2  &&  priorityHad1 == priorityHad2  &&  UseRandom::rndbool() ) {
		priorityHad2 = 0;
	}

	Lorentz5Momentum pClu = ptr->momentum();
	bool secondHad = false;
	Axis uSmear_v3;
	if ( priorityHad1  ||  priorityHad2 ) {

		double cluSmear;
		Lorentz5Momentum pQ;
		Lorentz5Momentum pQbar;
		if ( priorityHad1 > priorityHad2 ) {
			pQ    = ptr1->momentum();
			pQbar = ptr2->momentum();
			cluSmear = cluSmearHad1;
		} else {
			pQ    = ptr2->momentum();
			pQbar = ptr1->momentum();
			cluSmear = cluSmearHad2;
			secondHad = true;
		}

		// To find the momenta of the two hadrons in the parent cluster frame
		// we proceed as follows. First, we determine the unit vector parallel
		// to the direction of the constituent in the cluster frame. Then we
		// have to smear such direction using the following prescription:
		//  --- in theta angle w.r.t. such direction (not the z-axis),
		//      the drawing of the cosine of such angle is done via:
		//                 1.0 + cluSmear*log( rnd() )
		//      (repeated until it gives a value above -1.0)
		//  --- in phi angle w.r.t. such direction, the drawing is simply flat.
		// Then, given the direction in the parent cluster frame of one of the
		// two hadrons, it is straighforward to get the momenta of both hadrons
		// (always in the same parent cluster frame).

		pQ.boost( -pClu.boostVector() );    // boost from Lab to Cluster frame
		uSmear_v3 = pQ.vect().unit(); // Direction of the priority quark in COM frame
		// skip if cluSmear is too small
		if ( cluSmear > 0.001 ) {
			// generate the smearing angle
			double cosSmear;
			do cosSmear = 1.0 + cluSmear*log( UseRandom::rnd() );
			while ( cosSmear < -1.0 );
			double sinSmear = sqrt( 1.0 - sqr(cosSmear) );
			// calculate rotation to z axis
			LorentzRotation rot;
			double sinth(sqrt(1.-sqr(uSmear_v3.z())));
			if(abs(uSmear_v3.perp2()/uSmear_v3.z())>1e-10)
				rot.setRotate(-acos(uSmear_v3.z()),
						Axis(-uSmear_v3.y()/sinth,uSmear_v3.x()/sinth,0.));
			// + random azimuthal rotation
			rot.rotateZ(UseRandom::rnd()*twopi);
			// set direction in rotated frame
			Lorentz5Vector<double> ltemp(0.,sinSmear,cosSmear,0.);
			// rotate back
			rot.invert();
			ltemp *= rot;
			uSmear_v3 = ltemp.vect();
		}
	}
	else {

		// Isotropic decay: flat in cosTheta and phi.
		uSmear_v3 = Axis(1.0, 0.0, 0.0);  // just to set the rho to 1
		uSmear_v3.setTheta( acos( UseRandom::rnd( -1.0 , 1.0 ) ) );
		uSmear_v3.setPhi( UseRandom::rnd( -pi , pi ) );

	}

	pair<tcPDPtr,tcPDPtr> dataPair
		= _hadronSpectrum->chooseHadronPair(ptr->mass(),
				ptr1data,
				ptr2data);
	if(dataPair.first  == tcPDPtr() ||
			dataPair.second == tcPDPtr()) return pair<PPtr,PPtr>();

	// Create the two hadron particle objects with the specified id.
	PPtr ptrHad1,ptrHad2;
	// produce the hadrons on mass shell
	if(_onshell) {
		ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
		ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
	}
	// produce the hadrons with mass given by the mass generator
	else {
		unsigned int ntry(0);
		do {
			ptrHad1 = dataPair.first ->produceParticle();
			ptrHad2 = dataPair.second->produceParticle();
			++ntry;
		}
		while(ntry<_masstry&&ptrHad1->mass()+ptrHad2->mass()>ptr->mass());
		// if fails produce on shell and issue warning (should never happen??)
		if( ptrHad1->mass() + ptrHad2->mass() > ptr->mass() ) {
			generator()->log() << "Failed to produce off-shell hadrons in "
				<< "ClusterDecayer::decayIntoTwoHadrons producing hadrons "
				<< "on-shell" << endl;
			ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
			ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
		}
	}

	if (!ptrHad1  ||  !ptrHad2) {
		ostringstream s;
		s << "ClusterDecayer::decayIntoTwoHadrons ***Cannot create the two hadrons***"
			<< dataPair.first ->PDGName() << " and "
			<< dataPair.second->PDGName();
		cerr << s.str() << endl;
		generator()->logWarning( Exception(s.str(), Exception::warning) );
	} else {

		Lorentz5Momentum pHad1, pHad2;  // 5-momentum vectors that we have to determine
		if ( secondHad ) uSmear_v3 *= -1.0;

		if (pClu.m() < ptrHad1->mass()+ptrHad2->mass() ) {
			throw Exception() << "Impossible Kinematics in ClusterDecayer::decayIntoTwoHadrons()"
				<< Exception::eventerror;
		}
		// 5-momentum vectors that we have to determine
		Kinematics::twoBodyDecay(pClu,ptrHad1->mass(),ptrHad2->mass(),uSmear_v3,
				pHad1,pHad2);
		ptrHad1->set5Momentum(pHad1);
		ptrHad2->set5Momentum(pHad2);

		// Determine the positions of the two children clusters.
		LorentzPoint positionHad1 = LorentzPoint();
		LorentzPoint positionHad2 = LorentzPoint();
		calculatePositions(pClu, ptr->vertex(), pHad1, pHad2, positionHad1, positionHad2);
		ptrHad1->setVertex(positionHad1);
		ptrHad2->setVertex(positionHad2);
	}
	return pair<PPtr,PPtr>(ptrHad1,ptrHad2);
}

pair<PPtr,PPtr> ClusterDecayer::decayIntoTwoHadronsNew(tClusterPtr ptr) {
  using Constants::pi;
  using Constants::twopi;
	/* New test for cluster decay modelling according to matrix element
	 * C-> q,qbar->h1,h2  using |M|^2=1/((p1-p3)^2-(m1-m3)^2)^2 which is Tchannel-like
	 * */
  if ( ptr->numComponents() != 2 ) {
    generator()->logWarning( Exception("ClusterDecayer::decayIntoTwoHadronsNew "
				       "***Still cluster with not exactly 2 components*** ",
				       Exception::warning) );
    return pair<PPtr,PPtr>();
  }

  // Extract the id and particle pointer of the two components of the cluster.
  tPPtr ptr1 = ptr->particle(0);
  tPPtr ptr2 = ptr->particle(1);
  tcPDPtr ptr1data = ptr1->dataPtr();
  tcPDPtr ptr2data = ptr2->dataPtr();

  Lorentz5Momentum pClu = ptr->momentum();

  pair<tcPDPtr,tcPDPtr> dataPair
    = _hadronSpectrum->chooseHadronPair(ptr->mass(),
					 ptr1data,
					 ptr2data);
  if(dataPair.first  == tcPDPtr() ||
     dataPair.second == tcPDPtr()) return pair<PPtr,PPtr>();

  // Create the two hadron particle objects with the specified id.
  PPtr ptrHad1,ptrHad2;
  // produce the hadrons on mass shell
	if(_onshell) {
		ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
		ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
	}
	// produce the hadrons with mass given by the mass generator
	else {
		unsigned int ntry(0);
		do {
			ptrHad1 = dataPair.first ->produceParticle();
			ptrHad2 = dataPair.second->produceParticle();
			++ntry;
		}
		while(ntry<_masstry&&ptrHad1->mass()+ptrHad2->mass()>ptr->mass());
		// if fails produce on shell and issue warning (should never happen??)
		if( ptrHad1->mass() + ptrHad2->mass() > ptr->mass() ) {
			generator()->log() << "Failed to produce off-shell hadrons in "
				<< "ClusterDecayer::decayIntoTwoHadronsNew producing hadrons "
				<< "on-shell" << endl;
			ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
			ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
		}
	}

  if (!ptrHad1  ||  !ptrHad2) {
    ostringstream s;
    s << "ClusterDecayer::decayIntoTwoHadronsNew ***Cannot create the two hadrons***"
      << dataPair.first ->PDGName() << " and "
      << dataPair.second->PDGName();
    cerr << s.str() << endl;
    generator()->logWarning( Exception(s.str(), Exception::warning) );
  } else {

	  if (pClu.m() < ptrHad1->mass()+ptrHad2->mass() ) {
		  throw Exception() << "Impossible Kinematics in ClusterDecayer::decayIntoTwoHadronsNew()"
			  << Exception::eventerror;
	  }
	  Lorentz5Momentum pHad1,pHad2;
	  Axis uSmear_v3;
		Lorentz5Momentum pQ = ptr1->momentum();

		pQ.boost( -pClu.boostVector() );    // boost from Lab to Cluster frame
		const Axis dirQ = pQ.vect().unit(); // Direction of the quark in COM frame

	  switch (abs(_kinematics))
	  {
		  case 1:
			  {
				  Energy mHad1   = ptrHad1->mass();
				  Energy mHad2   = ptrHad2->mass();
				  Energy mConst1 = ptr1data->constituentMass();
				  Energy mConst2 = ptr2data->constituentMass();
					// New sampling from distribution 1/((p1-p3)^2-MS^2)^2
					// with MS = max((mHad1-mConst1),(mHad2-mConst2))
					// Maximum of masses 
					Energy mConst = mConst2;
					Energy mHad   = mHad2;
					if (sqr(mHad1-mConst1) > sqr(mHad2-mConst2)) {
						mConst = mConst1;
						mHad   = mHad1;
					}
					double PcomIniDivM = Kinematics::pstarTwoBodyDecay(pClu.m(),mConst1,mConst2)/mConst;
					double PcomFinDivM = Kinematics::pstarTwoBodyDecay(pClu.m(),mHad1,mHad2)/mHad;
					// A = (EcomIni1*EcomFin1-mHad*mConst)/PcomIni*PcomFin;
					// Note expm1 and log1p for more numerically stable results
					double Ainv = (PcomIniDivM*PcomFinDivM)/(expm1(0.5*log1p(PcomIniDivM*PcomIniDivM + PcomFinDivM*PcomFinDivM + PcomFinDivM*PcomFinDivM*PcomIniDivM*PcomIniDivM )));
					double cosTheta,phi;
					do {
						uSmear_v3 = dirQ;
						cosTheta = Kinematics::sampleCosTchannel(Ainv);
						// If no change in angle keep the direction fixed
						if (fabs(cosTheta-1.0)>1e-14) {
							// rotate to sampled angles
							uSmear_v3.rotate(acos(cosTheta),dirQ.orthogonal());
							phi  = UseRandom::rnd(-pi,pi);
							uSmear_v3.rotate(phi,dirQ);
						}
					}
					// filter out not corresponding hemisphere
					while ( _kinematics<0 && dirQ.dot(uSmear_v3)<=0.0);
					break;
			  }
		  case 2:
			  {
					// uSmear_v3 changed to isotropic direction
					Axis iso;
					do {
						iso = Axis(1.0, 0.0, 0.0);  // just to set the rho to 1
						iso.setTheta( acos( UseRandom::rnd( -1.0 , 1.0 ) ) );
						iso.setPhi( UseRandom::rnd( -pi , pi ) );
					}
					// filter out not corresponding hemisphere
					while ( _kinematics<0 && dirQ.dot(iso)<=0.0);
					uSmear_v3 = iso;
				  break;
			  }
		  case 3:
			  {
					// set to be aligned
					uSmear_v3 = dirQ;
				  break;
			  }
		  default:
			  assert(false);
	  }
		// sanity check
		if (_kinematics<0 && !(uSmear_v3.dot(dirQ)>0.0)) {
			std::cout << "uSmear_v3.dirQ  "<< uSmear_v3.dot(dirQ) << std::endl;
		}
		assert(_kinematics>0 || (_kinematics<0 && uSmear_v3.dot(dirQ)>0.0));
		// 5-momentum vectors that we have to determine
		Kinematics::twoBodyDecay(pClu,ptrHad1->mass(),ptrHad2->mass(),uSmear_v3,
				pHad1,pHad2);
	  ptrHad1->set5Momentum(pHad1);
	  ptrHad2->set5Momentum(pHad2);

	  // Determine the positions of the two children clusters.
	  LorentzPoint positionHad1 = LorentzPoint();
	  LorentzPoint positionHad2 = LorentzPoint();
	  calculatePositions(pClu, ptr->vertex(), pHad1, pHad2, positionHad1, positionHad2);
	  ptrHad1->setVertex(positionHad1);
	  ptrHad2->setVertex(positionHad2);
  }
  return pair<PPtr,PPtr>(ptrHad1,ptrHad2);
}


void ClusterDecayer::
calculatePositions(const Lorentz5Momentum &pClu,
		   const LorentzPoint &positionClu,
		   const Lorentz5Momentum &,
		   const Lorentz5Momentum &,
		   LorentzPoint &positionHad1,
		   LorentzPoint &positionHad2 ) const {
  // First, determine the relative positions of the children hadrons
  // with respect to their parent cluster, in the cluster reference frame,
  // assuming gaussian smearing with width inversely proportional to the
  // parent cluster mass.
  Length smearingWidth = hbarc / pClu.m();
  LorentzDistance distanceHad[2];
  for ( int i = 0; i < 2; i++ ) {   // i==0 is the first hadron; i==1 is the second one
    Length delta[4]={ZERO,ZERO,ZERO,ZERO};
    // smearing of the four components of the LorentzDistance, two at the same time to improve speed
    for ( int j = 0; j < 3; j += 2 ) {
      delta[j] = UseRandom::rndGauss(smearingWidth, Length(ZERO));
      delta[j+1] = UseRandom::rndGauss(smearingWidth, Length(ZERO));
    }
    // set the distance
    delta[0] = abs(delta[0]) +sqrt(sqr(delta[1])+sqr(delta[2])+sqr(delta[3]));
    distanceHad[i] = LorentzDistance(delta[1],delta[2],delta[3],delta[0]);
    // Boost such relative positions of the children hadrons,
    // with respect to their parent cluster,
    // from the cluster reference frame to the Lab frame.
    distanceHad[i].boost(pClu.boostVector());
  }
  // Finally, determine the absolute positions of the children hadrons
  // in the Lab frame.
  positionHad1 = distanceHad[0] + positionClu;
  positionHad2 = distanceHad[1] + positionClu;
}
