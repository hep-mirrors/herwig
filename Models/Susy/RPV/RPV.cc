// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPV class.
//

#include "RPV.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RPV::persistentOutput(PersistentOStream & os ) const {
  os << lambdaLLE_ << lambdaLQD_ << lambdaUDD_ << ounit(vnu_,GeV)
     << upSquarkMix_ << downSquarkMix_ << triLinearOnly_
     << LLEVertex_ << LQDVertex_ << UDDVertex_ 
     << ounit(epsilon_,GeV) << ounit(epsB_,GeV);
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_ >> iunit(vnu_,GeV)
     >> upSquarkMix_ >> downSquarkMix_ >> triLinearOnly_
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_ 
     >> iunit(epsilon_,GeV) >> iunit(epsB_,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<RPV,MSSM>
describeHerwigRPV("Herwig::RPV", "HwSusy.so HwRPV.so");

void RPV::Init() {

  static ClassDocumentation<RPV> documentation
    ("The RPV class is the base class for the implementation of the"
     " R-parity violating MSSM.");

  static Reference<RPV,AbstractFFSVertex> interfaceLLEVertex
    ("Vertex/LLE",
     "The vertex for the trillinear LLE interaction",
     &RPV::LLEVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceLQDVertex
    ("Vertex/LQD",
     "The vertex for the trillinear LQD interaction",
     &RPV::LQDVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceUDDVertex
    ("Vertex/UDD",
     "The vertex for the trillinear UDD interaction",
     &RPV::UDDVertex_, false, false, true, false, false);

  static Switch<RPV,bool> interfaceTriLinearOnly
    ("TriLinearOnly",
     "Only include trilinears and take rest of model to be MSSM",
     &RPV::triLinearOnly_, false, false, false);
  static SwitchOption interfaceTriLinearOnlyYes
    (interfaceTriLinearOnly,
     "Yes",
     "Trilinears + MSSM",
     true);
  static SwitchOption interfaceTriLinearOnlyNo
    (interfaceTriLinearOnly,
     "No",
     "All RPV couplings and mixings",
     false);

}

void RPV::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm != 0) 
      throw Exception() << "R-parity violating MSSM model"
			<< " used but NMSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 1) throw Exception() << "RPV model used but no RPV in input file"
				    << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "RPV model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "RPV model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the RPV parameters
  // lambda
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("rvlamlle");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLLE_[i][j][k] = it->second;
    }
  }
  // lambda'
  pit=parameters().find("rvlamlqd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLQD_[i][j][k] = it->second;
    }
  }
  // lambda''
  pit=parameters().find("rvlamudd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaUDD_[i][j][k] = it->second;
    }
  }
  // sneutrino vevs
  pit=parameters().find("rvsnvev");
  vnu_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	vnu_[it->first-1] = it->second*GeV;
      }
    }
  }
  // bilinears
  pit=parameters().find("rvkappa");
  epsilon_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	epsilon_[it->first-1] = it->second*GeV;
      }
    }
  }
  // blinear soft terms
  pit=parameters().find("rvd");
  epsB_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	epsB_[it->first-1] = it->second*GeV;
      }
    }
  }
}

void RPV::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "rvamix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      CPoddHiggsMix(temp);
    }
    else if (name == "rvlmix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      ChargedHiggsMix(temp);
    }
    else if (name == "dsqmix" ) {
      createMixingMatrix(downSquarkMix_,name,it->second.second,it->second.first);
    }
    else if (name == "usqmix" ) {
      createMixingMatrix(upSquarkMix_,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
  // now adjust the mixing matrices to have our structure
  // first bloodly SPHENO as it doesn't obey the SLHA  
  map<string,StringMap>::const_iterator sit = info().find("spinfo");
  string program;
  if(sit!=info().end()) {
    StringMap::const_iterator pit = sit->second.find(1);
    if(pit!=sit->second.end()) program = pit->second;
  }
  if(program=="SPheno") {
    map<string,ParamMap>::const_iterator fit=parameters().find("mass");
    if(fit==parameters().end()) 
      throw Exception() << "BLOCK MASS not found in input file"
			<< " can't set masses of SUSY particles"
			<< Exception::runerror;
    // adjust the charged scalars
    map<double,long> massMap;
    massMap[findValue(fit,     37,"mass",     "37")] =      37;
    massMap[findValue(fit,1000011,"mass","1000011")] = 1000011;
    massMap[findValue(fit,1000013,"mass","1000013")] = 1000013;
    massMap[findValue(fit,1000015,"mass","1000015")] = 1000015;
    massMap[findValue(fit,2000011,"mass","2000011")] = 2000011;
    massMap[findValue(fit,2000013,"mass","2000013")] = 2000013;
    massMap[findValue(fit,2000015,"mass","2000015")] = 2000015;
    massMap[getParticleData(24)->mass()/GeV        ] =      24;
    vector<int> move;
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     37) move.push_back(0);
      else if(mit->second==1000011) move.push_back(1);
      else if(mit->second==1000013) move.push_back(2);
      else if(mit->second==1000015) move.push_back(3);
      else if(mit->second==2000011) move.push_back(4);
      else if(mit->second==2000013) move.push_back(5);
      else if(mit->second==2000015) move.push_back(6);
      else if(mit->second==     24) move.push_back(7);
    }
    CMatrix oldMat = ChargedHiggsMix()->getMatrix();
    CMatrix newMat(8,vector<Complex>(8,0.));
    for(unsigned int ix=0;ix<8;++ix) {
      for(unsigned int iy=0;iy<8;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    ChargedHiggsMix(new_ptr(MixingMatrix(newMat,ChargedHiggsMix()->getIds())));
    // adjust the pseudoscalars
    massMap.clear();
    massMap[findValue(fit,     36,"mass",     "36")] =      36;
    // extract the pseudoscalar masses and change the ids for the pseudoscalars
    // to those from the SLHA if needed
    if(fit->second.find(2000012)!=fit->second.end()) {
      massMap[findValue(fit,2000012,"mass","2000012")] = 1000017;
      idMap().insert(make_pair(2000012,1000017));
    }
    else {
      massMap[findValue(fit,1000017,"mass","1000017")] = 1000017;
    }
    if(fit->second.find(2000014)!=fit->second.end()) {
      massMap[findValue(fit,2000014,"mass","2000014")] = 1000018;
      idMap().insert(make_pair(2000014,1000018));
    }
    else {
      massMap[findValue(fit,1000018,"mass","1000018")] = 1000018;
    }
    if(fit->second.find(2000016)!=fit->second.end()) {
      massMap[findValue(fit,2000016,"mass","2000016")] = 1000019;
      idMap().insert(make_pair(2000016,1000019));
    }
    else {
      massMap[findValue(fit,1000019,"mass","1000019")] = 1000019;
    }
    move.clear(); move.push_back(4);
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     36) move.push_back(0);
      else if(mit->second==1000017) move.push_back(1);
      else if(mit->second==1000018) move.push_back(2);
      else if(mit->second==1000019) move.push_back(3);
    }
    oldMat = CPoddHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPoddHiggsMix(new_ptr(MixingMatrix(newMat,CPoddHiggsMix()->getIds())));
    // adjust the neutral scalars
    massMap.clear();
    massMap[findValue(fit,     25,"mass",     "25")] =      25;
    massMap[findValue(fit,     35,"mass",     "35")] =      35;
    massMap[findValue(fit,1000012,"mass","1000012")] = 1000012;
    massMap[findValue(fit,1000014,"mass","1000014")] = 1000014;
    massMap[findValue(fit,1000016,"mass","1000016")] = 1000016;
    move.clear();
    for(map<double,long>::iterator mit=massMap.begin();mit!=massMap.end();++mit) {
      if     (mit->second==     25) move.push_back(0);
      else if(mit->second==     35) move.push_back(1);
      else if(mit->second==1000012) move.push_back(2);
      else if(mit->second==1000014) move.push_back(3);
      else if(mit->second==1000016) move.push_back(4);
    }
    oldMat = CPevenHiggsMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][iy] = oldMat[ix][iy];
    }
    CPevenHiggsMix(new_ptr(MixingMatrix(newMat,CPevenHiggsMix()->getIds())));
    // neutralino mixing
    move.resize(7);
    move[0] = 3; move[1] = 4; move[2] = 5; move[3] = 6;
    move[4] = 0; move[5] = 1; move[6] = 2;
    oldMat = neutralinoMix()->getMatrix();
    newMat = CMatrix(7,vector<Complex>(7,0.));
    for(unsigned int ix=0;ix<7;++ix) {
      for(unsigned int iy=0;iy<7;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
    // chargino mixing
    move.resize(5);
    move[0] = 3; move[1] = 4;
    move[2] = 0; move[3] = 1; move[4] = 2;
    oldMat = charginoUMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
    oldMat = charginoVMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[ix][move[iy]] = oldMat[ix][iy];
    }
    charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
  }
  // we don't want neutrinos first then neutralinos so swap them
  // neutralinos first then neutrinos
  if ( neutralinoMix()->size().first == 7 ) {
    vector<int> move(7);
    move[0] = 4; move[1] = 5; move[2] = 6;
    move[3] = 0; move[4] = 1; move[5] = 2; move[6] = 3;
    CMatrix oldMat = neutralinoMix()->getMatrix();
    CMatrix newMat(7,vector<Complex>(7,0.));
    for(unsigned int ix=0;ix<7;++ix) {
      for(unsigned int iy=0;iy<7;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    neutralinoMix(new_ptr(MixingMatrix(newMat,neutralinoMix()->getIds())));
  }
  // charginos the same, i.e. charginos first then charged leptons
  if(charginoUMix()->size().first  != charginoVMix()->size().first || 
     charginoUMix()->size().second != charginoVMix()->size().second )
    throw Exception() << "Chargino U and V mixing matrices must have the same size.\n"
		      << "Check your SLHA file!" << Exception::runerror;
  if ( charginoUMix()->size().first == 5 ) {
    vector<int> move(5);
    move[0] = 2; move[1] = 3; move[2] = 4;
    move[3] = 0; move[4] = 1;
    CMatrix oldMat = charginoUMix()->getMatrix();
    CMatrix newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    charginoUMix(new_ptr(MixingMatrix(newMat,charginoUMix()->getIds())));
    oldMat = charginoVMix()->getMatrix();
    newMat = CMatrix(5,vector<Complex>(5,0.));
    for(unsigned int ix=0;ix<5;++ix) {
      for(unsigned int iy=0;iy<5;++iy)
	newMat[move[ix]][move[iy]] = oldMat[ix][iy];
    }
    charginoVMix(new_ptr(MixingMatrix(newMat,charginoVMix()->getIds())));
  }

  const MatrixSize & n = neutralinoMix()->size();
  const MatrixSize & u = charginoUMix()->size();
  const MatrixSize & h = CPevenHiggsMix()->size();
  const MatrixSize & a = CPoddHiggsMix()->size();
  const MatrixSize & l = ChargedHiggsMix()->size();

  bool nBig   = n.first == 7 && n.second == 7;
  bool nSmall = n.first == 4 && n.second == 4;
  if ( ! (nBig || nSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Nmix " << n.first << ',' << n.second << '\n'
      << Exception::runerror; 

  bool uBig   = u.first == 5 && u.second == 5;
  bool uSmall = u.first == 2 && u.second == 2;
  if ( ! (uBig || uSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Umix " << u.first << ',' << u.second << '\n'
      << Exception::runerror; 

  bool hBig   = h.first == 5 && h.second == 5;
  bool hSmall = h.first == 2 && h.second == 2;
  if ( ! (hBig || hSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Hmix " << h.first << ',' << h.second << '\n'
      << Exception::runerror; 

  bool aBig = (a.first == 4 || a.first == 5) && a.second == 5;
  bool aSmall = a.first == 1 && a.second == 2;
  if ( ! (aBig || aSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "RVAmix " << a.first << ',' << a.second << '\n'
      << Exception::runerror; 

  bool lBig = (l.first == 7 || l.first == 8) && l.second == 8;
  bool lSmall = l.first == 1 && l.second == 2;
  if ( ! (lBig || lSmall) )
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "RVLmix " << l.first << ',' << l.second << '\n'
      << Exception::runerror; 


  bool allBig   = nBig && uBig && hBig && aBig && lBig;
  bool allSmall = nSmall && uSmall && hSmall && aSmall && lSmall;

  bool allSmallExceptN = nBig && uSmall && hSmall && aSmall && lSmall;
  
  if ( allSmallExceptN ) {
    cerr << "Warning: Truncating Nmix to 4,4 for consistency "
	 << "with other mixing matrices.\n";
    CMatrix oldMat = neutralinoMix()->getMatrix();
    CMatrix newMat(4,vector<Complex>(4,0.));
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	newMat[ix][iy] = oldMat[ix][iy];
    assert( neutralinoMix()->getIds().size() >= 4 );
    vector<long>::const_iterator beg = neutralinoMix()->getIds().begin();
    vector<long>::const_iterator end = beg + 4;
    neutralinoMix(new_ptr(MixingMatrix(newMat,vector<long>(beg,end))));
    return;
  }
  else if ( ! (allBig || allSmall) ) {
    throw Exception() 
      << "Mixing matrices have inconsistent sizes:\n"	
      << "(RV)Nmix " << n.first << ',' << n.second << '\n'
      << "(RV)Umix " << u.first << ',' << u.second << '\n'
      << "(RV)Hmix " << h.first << ',' << h.second << '\n'
      << "RVAmix   " << a.first << ',' << a.second << '\n'
      << "RVLmix   " << l.first << ',' << l.second << '\n'
      << Exception::runerror;
  }
  // reduce to MSSM + trilinear if requested
  if(triLinearOnly_) {
    // reduce size of n
    if(neutralinoMix()->size().first ==7) {
      CMatrix mix(4,vector<Complex>(4,0.));
      unsigned int irow=0;
      int imax[7]={-1,-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<7;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<7;++iy) {
       	  double value = abs((*neutralinoMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// neutralino
	if(imax[ix]<=3) {
	  for(unsigned int iy=0;iy<4;++iy) mix[irow][iy] = (*neutralinoMix())(ix,iy);
	  ++irow;
	  assert(irow<=4);
	}
	// neutrino
	else {
	  idMap()[neutralinoMix()->getIds()[ix]] = neutralinoMix()->getIds()[imax[ix]];
	}
      }
      vector<long> ids = neutralinoMix()->getIds();
      ids.resize(4);
      neutralinoMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of u
    if(charginoUMix()->size().first ==5) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*charginoUMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// chargino
       	if(imax[ix]<=1) {
     	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*charginoUMix())(ix,iy);
     	  ++irow;
     	  assert(irow<=2);
      	}
	// charged lepton
      	else {
      	  idMap()[abs(charginoUMix()->getIds()[ix])] = abs(charginoUMix()->getIds()[imax[ix]]);
     	}
      }
      vector<long> ids = charginoUMix()->getIds();
      ids.resize(2);
      charginoUMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of v
    if(charginoVMix()->size().first ==5) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*charginoVMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
	// chargino
       	if(imax[ix]<=1) {
     	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*charginoVMix())(ix,iy);
     	  ++irow;
     	  assert(irow<=2);
      	}
      }
      vector<long> ids = charginoVMix()->getIds();
      ids.resize(2);
      charginoVMix(new_ptr(MixingMatrix(mix,ids)));
    }
    // reduce size of pseudo scalar mixing
    if(CPoddHiggsMix()) {
      MixingVector hmix;
      double beta = atan(tanBeta());
      hmix.push_back(MixingElement(1,1,sin(beta)));
      hmix.push_back(MixingElement(1,2,cos(beta)));
      vector<long> ids(1,36);
      MixingMatrixPtr newMix = new_ptr(MixingMatrix(1,2));
      (*newMix).setIds(ids);
      for(unsigned int ix=0; ix < hmix.size(); ++ix)
	(*newMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
      CPoddHiggsMix(newMix);
    }
    // reduce size of true scalar mixing
    if(CPevenHiggsMix()->size().first==5) {
      double alpha(0.);
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int irow=0;
      int imax[5]={-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<5;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<5;++iy) {
       	  double value = abs((*CPevenHiggsMix())(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
      	}
       	// neutral Higgs
       	if(imax[ix]<=1) {
      	  for(unsigned int iy=0;iy<2;++iy) mix[irow][iy] = (*CPevenHiggsMix())(ix,iy);
	  if(irow==0)
	    alpha = atan2(-(*CPevenHiggsMix())(ix,0).real(),(*CPevenHiggsMix())(ix,1).real());
      	  ++irow;
      	  assert(irow<=2);
       	}
      }
      vector<long> ids = CPevenHiggsMix()->getIds();
      ids.resize(2);
      CPevenHiggsMix(new_ptr(MixingMatrix(mix,ids)));
      higgsMixingAngle(alpha);
    }
    else {
      bool readAlpha = false;
      map<string,ParamMap>::const_iterator pit=parameters().find("alpha");
      if(pit!=parameters().end()) {
	ParamMap::const_iterator it = pit->second.find(1);
	if(it!=pit->second.end()) {
	  readAlpha = true;
	  higgsMixingAngle(it->second);
	}
      }
      if(!readAlpha) 
	throw Exception() << "In the RPV model BLOCK ALPHA which must be"
			  << " present in the SLHA file is missing"
			  << Exception::runerror;
    }
    // reduce size of charged scalar mixing
    if(ChargedHiggsMix()->size().first>=7) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int istau(0);
      int imax[8]={-1,-1,-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<ChargedHiggsMix()->size().first;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<8;++iy) {
	  double value = abs((*ChargedHiggsMix())(ix,iy));
	  if(value>maxComp) {
	    maxComp = value;
	    imax[ix] = iy;
	  }
	}
	if(imax[ix]<=1) imax[ix]=0;
	else --imax[ix];
	idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = abs(ChargedHiggsMix()->getIds()[imax[ix]]);
      	if(abs(ChargedHiggsMix()->getIds()[imax[ix]])%10==5) {
      	  mix[istau][0] = (*ChargedHiggsMix())(ix,4);
      	  mix[istau][1] = (*ChargedHiggsMix())(ix,7);
	  if(istau==0) idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = 1000015;
	  else         idMap()[abs(ChargedHiggsMix()->getIds()[ix])] = 2000015;
      	  istau+=1;
      	  assert(istau<=2);
      	}
      }
      // set up the stau mixing matrix
      vector<long> ids(2);
      ids[0] = 1000015;
      ids[1] = 2000015;
      stauMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 7x7
      MixingVector hmix;
      double beta = atan(tanBeta());
      hmix.push_back(MixingElement(1,1,sin(beta)));
      hmix.push_back(MixingElement(1,2,cos(beta)));
      ids.clear();
      ids.resize(1,37);
      MixingMatrixPtr newMix = new_ptr(MixingMatrix(1,2));
      (*newMix).setIds(ids);
      for(unsigned int ix=0; ix < hmix.size(); ++ix)
      	(*newMix)(hmix[ix].row-1,hmix[ix].col-1) = hmix[ix].value;
      ChargedHiggsMix(newMix);
    }
    // reduce size of up squark mixing
    if( upSquarkMix_ ) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int istop(0);
      int imax[6]={-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<6;++ix) {
	double maxComp(0.);
	for(unsigned int iy=0;iy<6;++iy) {
	  double value = abs((*upSquarkMix_)(ix,iy));
	  if(value>maxComp) {
	    maxComp = value;
	    imax[ix] = iy;
	  }
	}
	idMap()[upSquarkMix_->getIds()[ix]] = upSquarkMix_->getIds()[imax[ix]];
	if(upSquarkMix_->getIds()[imax[ix]]%10==6) {
	  mix[istop][0] = (*upSquarkMix_)(ix,2);
	  mix[istop][1] = (*upSquarkMix_)(ix,5);
	  if(istop==0) idMap()[upSquarkMix_->getIds()[ix]] = 1000006;
	  else         idMap()[upSquarkMix_->getIds()[ix]] = 2000006;
	  istop+=1;
	  assert(istop<=2);
	}
      }
      // set up the stop mixing matrix
      vector<long> ids(2);
      ids[0] = 1000006;
      ids[1] = 2000006;
      stopMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 6x6
      upSquarkMix_ = MixingMatrixPtr();
    }
    // reduce size of down squark mixing
    if( downSquarkMix_ ) {
      CMatrix mix(2,vector<Complex>(2,0.));
      unsigned int isbot(0);
      int imax[6]={-1,-1,-1,-1,-1,-1};
      for(unsigned int ix=0;ix<6;++ix) {
      	double maxComp(0.);
      	for(unsigned int iy=0;iy<6;++iy) {
      	  double value = abs((*downSquarkMix_)(ix,iy));
       	  if(value>maxComp) {
       	    maxComp = value;
       	    imax[ix] = iy;
       	  }
       	}
      	idMap()[downSquarkMix_->getIds()[ix]] = downSquarkMix_->getIds()[imax[ix]];
      	if(downSquarkMix_->getIds()[imax[ix]]%10==5) {
      	  mix[isbot][0] = (*downSquarkMix_)(ix,2);
      	  mix[isbot][1] = (*downSquarkMix_)(ix,5);
	  if(isbot==0) idMap()[downSquarkMix_->getIds()[ix]] = 1000005;
	  else         idMap()[downSquarkMix_->getIds()[ix]] = 2000005;
      	  isbot+=1;
      	  assert(isbot<=2);
      	}
      }
      // set up the sbottom mixing matrix
      vector<long> ids(2);
      ids[0] = 1000005;
      ids[1] = 2000005;
      sbottomMix(new_ptr(MixingMatrix(mix,ids)));
      // delete 6x6
      downSquarkMix_ = MixingMatrixPtr();
    }
    // get the masses as we will have to mess with them
    map<string,ParamMap>::const_iterator fit=parameters().find("mass");
    ParamMap theMasses = fit->second;
    for(long id=1000017;id<=1000019;++id) {
      if(theMasses.find(id)==theMasses.end()) continue;
      double mass = abs(theMasses.find(id)->second);
      // find scalar partner
      double mdiff=1e30;
      long new_id = 0;
      for(ParamMap::const_iterator it=theMasses.begin();it!=theMasses.end();++it) {
	double diff = abs(abs(it->second)-mass);
	if(diff<mdiff) {
	  mdiff = diff;
	  new_id = it->first;
	}
      }
      if(idMap().find(new_id)!=idMap().end()) {
	idMap()[id] = idMap()[new_id];
      }
      else {
	idMap()[id] = new_id;
      }
    }
  }
}

void RPV::doinit() {
  MSSM::doinit();
  addVertex(LLEVertex_);
  addVertex(LQDVertex_);
  addVertex(UDDVertex_);
}
