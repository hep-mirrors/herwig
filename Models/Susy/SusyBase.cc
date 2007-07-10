// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SusyBase class.
//

#include "SusyBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"

using namespace Herwig;

void SusyBase::doinit() throw(InitException) {
  addVertex(theZSFSFVertex);
  addVertex(thePSFSFVertex);
  addVertex(theWSFSFVertex);
  addVertex(theNFSFVertex);
  addVertex(theGFSFVertex);
  addVertex(theHSFSFVertex);
  addVertex(theCFSFVertex);
  addVertex(theGSFSFVertex);
  addVertex(theGGSQSQVertex);
  addVertex(theGSGSGVertex);
  addVertex(theNNZVertex);
  addVertex(theCCPVertex);
  addVertex(theCCZVertex);
  addVertex(theCNWVertex);
  addVertex(theGOGOHVertex);
  addVertex(theWHHVertex);
  addVertex(theHHHVertex);
  StandardModel::doinit();
}

void SusyBase::persistentOutput(PersistentOStream & os) const {
  os << theNMix << theUMix << theVMix 
     << theZSFSFVertex << thePSFSFVertex << theWSFSFVertex 
     << theNFSFVertex << theGFSFVertex << theHSFSFVertex << theCFSFVertex 
     << theGSFSFVertex << theGGSQSQVertex 
     << theGSGSGVertex << theNNZVertex << theCCPVertex 
     << theCCZVertex << theCNWVertex << theSSFFHVertex << theGOGOHVertex
     << theSSWWHVertex << theWHHVertex << theHHHVertex << theSSHGGVertex
     << _tanbeta << ounit(_mu,GeV) 
     << ounit(theMone,GeV) << ounit(theMtwo,GeV) << ounit(theMthree,GeV);
}

void SusyBase::persistentInput(PersistentIStream & is, int) {
  is >> theNMix >> theUMix >> theVMix  
     >> theZSFSFVertex >> thePSFSFVertex >> theWSFSFVertex 
     >> theNFSFVertex >> theGFSFVertex >> theHSFSFVertex >> theCFSFVertex 
     >> theGSFSFVertex >> theGGSQSQVertex >> theGSGSGVertex 
     >> theNNZVertex >> theCCPVertex >> theCCZVertex >> theCNWVertex
     >> theSSFFHVertex >> theGOGOHVertex >> theSSWWHVertex >> theWHHVertex
     >> theHHHVertex >> theSSHGGVertex
     >> _tanbeta >> iunit(_mu,GeV) 
     >> iunit(theMone,GeV) >> iunit(theMtwo,GeV) >> iunit(theMthree,GeV);
}

ClassDescription<SusyBase> SusyBase::initSusyBase;
// Definition of the static class description member.

void SusyBase::Init() {

  static ClassDocumentation<SusyBase> documentation
    ("This is the base class for any SUSY model.");

  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexZSS
    ("Vertex/ZSFSF",
     "Reference to Susy ZSSVertex",
     &SusyBase::theZSFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexPSS
    ("Vertex/PSFSF",
     "Reference to Susy P SF SF vertex",
     &SusyBase::thePSFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexWSS
    ("Vertex/WSFSF",
     "Reference to Susy W SF SF vertex",
     &SusyBase::theWSFSFVertex, false, false, true, false);

  
  static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexNFSF
    ("Vertex/NFSF",
     "Reference to the neutralino-fermion-sfermion vertex",
     &SusyBase::theNFSFVertex, false, false, true, false);

  static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexGFSF
    ("Vertex/GFSF",
     "Reference to the gluino-fermion-sfermion vertex",
     &SusyBase::theGFSFVertex, false, false, true, false);
  
  static Reference<SusyBase,Helicity::SSSVertex> interfaceVertexHSFSF
    ("Vertex/HSFSF",
     "Reference to the Higgs-fermion-sfermion vertex",
     &SusyBase::theHSFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexCFSF
   ("Vertex/CFSF",
      "Reference to the chargino-fermion-sfermion vertex",
      &SusyBase::theCFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexGSFSF
   ("Vertex/GSFSF",
      "Reference to the gluon-sfermion-sfermion vertex",
      &SusyBase::theGSFSFVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VVSSVertex> interfaceVertexGGSS
     ("Vertex/GGSQSQ",
      "Reference to the gluon-gluon-squark-squark vertex",
      &SusyBase::theGGSQSQVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexGSGSG
     ("Vertex/GSGSG",
      "Reference to the gluon-gluino-gluino vertex",
      &SusyBase::theGSGSGVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexNNZ
    ("Vertex/NNZ",
     "Reference to Z-~chi_i0-~chi_i0 vertex",
     &SusyBase::theNNZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCCP
    ("Vertex/CCP",
     "Reference to ~chi_i+-~chi_i-photon vertex",
     &SusyBase::theCCPVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCCZ
    ("Vertex/CCZ",
     "Reference to ~chi_i+-~chi_i-Z vertex",
     &SusyBase::theCCZVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFVVertex> interfaceVertexCNW
    ("Vertex/CNW",
     "Reference to ~chi_i+-chi_i0-W vertex",
     &SusyBase::theCNWVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexSSFFH
   ("Vertex/SSFFH",
    "Reference to the fermion-antifermion-higgs vertex",
    &SusyBase::theSSFFHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::FFSVertex> interfaceVertexGOGOH
   ("Vertex/GOGOH",
    "Reference to the gaugino-gaugino-higgs vertex",
    &SusyBase::theGOGOHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VVSVertex> interfaceVertexWWH
   ("Vertex/SSWWH",
    "Reference to the boson-boson-higgs vertex",
    &SusyBase::theSSWWHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::VSSVertex> interfaceVertexWHH
    ("Vertex/SSWHH",
     "Reference to Susy WHHVertex",
     &SusyBase::theWHHVertex, false, false, true, false);

   static Reference<SusyBase,Helicity::SSSVertex> interfaceVertexHHH
    ("Vertex/HHH",
     "Triple higgs coupling",
     &SusyBase::theHHHVertex, false, false, true, false);
   
   static Reference<SusyBase,GeneralSVVVertex> interfaceVertexSSHGG
    ("Vertex/SSHGG",
     "The coupling of a higgs to 2 gluons",
     &SusyBase::theSSHGGVertex, false, false, true, false);

}

void SusyBase::readSetup(istream &is) throw(SetupException) {
  string name = dynamic_ptr_cast<istringstream*>(&is)->str();
  ifstream file(name.c_str());
  if(!file)
     throw SetupException() << "A problem occurred in opening the file: "
			    << name;
  string line;
  //function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  //Create a directory for the mixing matrices
  BaseRepository::CreateDirectory("/Defaults/MixingMatrices/");
  //Start to read in the file
  while(getline(file, line)) {
    // ignore comment lines
    if(line[0] == '#') continue;
    // make everything lower case
    transform(line.begin(), line.end(), line.begin(), pf);
    // start of a block
    if(line.find("block") == 0) {
      string name = StringUtils::car(StringUtils::cdr(line), " #");
      // mixing matrix
      if((name.find("mix")  != string::npos && 
	  name.find("hmix") != 0)||
	 name.find("au") == 0||name.find("ad") == 0||
	 name.find("ae") == 0 ) {
	unsigned int row(0),col(0);
	vector<MixingElement> vals = readMatrix(file,row,col);
	_mixings[name] = make_pair(make_pair(row,col),vals);
      }
      // decays
      else if( name.find("decay") == 0 ) {
	readDecay(file, line);
      }
      else if( name.find("info") == string::npos) {
	readBlock(file,name);
      }
    }
  }
  // extract the relevant parameters
  extractParameters();
  // create the mixing matrices we need
  createMixingMatrices();
  // set the masses, this has to be done after the 
  // mixing matrices
  resetRepositoryMasses();
}

void SusyBase::readBlock(ifstream & ifs,string name) throw(SetupException) {
  if(!ifs)
    throw SetupException() << "The input stream is in a bad state"
			   << Exception::runerror;
  string line;
  ParamMap store;
  // special for the alpha block
  if(name.find("alpha") == 0 ) {
    double alpha;
    getline(ifs, line);
    istringstream iss(line);
    iss >> alpha;
    store.insert(make_pair(1,alpha));
  }
  else {
    while(getline(ifs, line)) {
      if(line[0] == '#') {
	if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	    ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
	else continue;
      }
      istringstream is(line);
      long index;
      double value;   
      is >> index >> value;
      store.insert(make_pair(index, value));
      if(ifs.peek() == 'B' || ifs.peek() == 'D' ||
	 ifs.peek() == 'b' || ifs.peek() == 'd' ||
	 ifs.peek() == '#') break;
    }
  }
  _parameters[name]=store;
}

void SusyBase::readDecay(ifstream & ifs, 
			 string decay) const throw(SetupException){
  if(!ifs)
    throw SetupException() <<"The input stream is in a bad state";
  //if the particle is stable then next line will be the start of a 
  //new decaymode or a new block so just return if this is found
  if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
      ifs.peek() == 'd' || ifs.peek() == 'b' ) return;
  istringstream iss(decay);
  string dummy;
  long parent;
  double width;
  iss >> dummy >> parent >> width;
  //set the width
  tPDPtr inpart = getParticleData(parent);
  if(!inpart) throw SetupException() << "Decaying particle with PDG code "
				     << parent << " does not exist in "
				     << "SusyBase::readDecay()" 
				     << Exception::runerror;
  inpart->width(width*GeV);
  string tag = "decaymode " + inpart->PDGName() + "->";
  string line;
  while(getline(ifs, line)) {
    if(line[0] == '#') {
      if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	  ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
      else continue;
    }
    istringstream is(line);
    double brat(0.);
    unsigned int nda(0);
    is >> brat >> nda;
    vector<long> products(nda);
    vector<long>::size_type i(0);
    while(is && ++i < nda + 1)
      is >> products[i - 1];
    if( products.size() > 0 )
      createDecayMode(tag, products, brat);
     if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
         ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
  }
}

const vector<MixingElement>
SusyBase::readMatrix(ifstream & ifs, 
		     unsigned int & row, unsigned int & col) throw(SetupException) {
  if(!ifs)
    throw SetupException() << "The input stream is in a bad state";
  string line;
  unsigned int rowmax(0), colmax(0);
  vector<MixingElement> values;
  while(getline(ifs, line)) {
    if(line[0] == '#') {
      if( ifs.peek() == 'D' || ifs.peek() == 'B' ||
	  ifs.peek() == 'd' || ifs.peek() == 'b' ) break;
      else continue;
    }
    istringstream is(line);
    unsigned int index1, index2;
    double real(0.), imag(0.);   
    is >> index1 >> index2 >> real >> imag;
    values.push_back(MixingElement(index1,index2,Complex(real, imag)));
    if(index1 > rowmax) rowmax = index1; 
    if(index2 > colmax) colmax = index2; 
    if( ifs.peek() == 'B' || ifs.peek() == 'D' ||
        ifs.peek() == 'b' || ifs.peek() == 'd' ) break;
  }
  col=colmax;
  row=rowmax;
  return values;
}

void SusyBase::createDecayMode(string tag, vector<long> products,
			       double brat) const {
  ostringstream cmd;
  cmd << tag;
  vector<long>::size_type nda = products.size();
  for(vector<long>::size_type i = 0; i < nda; ++i) {
    tPDPtr prod = getParticleData(products[i]);
    if(!prod) throw SetupException() << "Unknown decay product with PDG code "
				     << products[i] << " in SusyBase::createDecayMode()"
				     << Exception::runerror;
    cmd << prod->PDGName();
    if(i != nda - 1) cmd << ",";
    else cmd << "; ";
  }
  cmd << brat << " 1 /Defaults/Decays/Mambo";
  Repository::exec(cmd.str(), cerr);
}

void SusyBase::createMixingMatrix(MixingMatrixPtr & matrix,
				  string name, const vector<MixingElement> & values,
				  pair<unsigned int,unsigned int> size) {
  string dir = "/Defaults/MixingMatrices/";
  MixingMatrixPtr mix = new_ptr(MixingMatrix(size.first,size.second));
  for(unsigned int ix=0; ix<values.size();++ix) {
    (*mix)(values[ix].row-1,values[ix].col-1) = values[ix].value;
  }
  if(name == "nmix") {
    matrix = mix;
    vector<long> ids(4);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035; 
    (*mix).setIds(ids);
  }
  else if(name == "nmnmix") {
    matrix = mix;
    vector<long> ids(5);
    ids[0] = 1000022; ids[1] = 1000023; 
    ids[2] = 1000025; ids[3] = 1000035;
    ids[4] = 1000045; 
    (*mix).setIds(ids);
  }
  else if(name == "umix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    (*mix).setIds(ids);
  }
  else if(name == "vmix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 1000024; ids[1] = 1000037; 
    (*mix).setIds(ids);
  }
  else if(name == "stopmix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 1000006; ids[1] = 2000006; 
    (*mix).setIds(ids);
  }
  else if(name == "sbotmix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 1000005; ids[1] = 2000005; 
    (*mix).setIds(ids);
  }
  else if(name == "staumix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 1000015; ids[1] = 2000015; 
    (*mix).setIds(ids);
  }
  else if(name == "nmhmix") {
    matrix = mix;
    vector<long> ids(3);
    ids[0] = 25; ids[1] = 35; ids[2] = 45;
    (*mix).setIds(ids);
  }
  else if(name == "nmamix") {
    matrix = mix;
    vector<long> ids(2);
    ids[0] = 36; ids[1] = 46;
    (*mix).setIds(ids);
  }
  else
    throw SetupException() << "Cannot find correct title for mixing matrix "
			   << name << Exception::runerror;
  BaseRepository::Register(mix, dir + name);
}

void SusyBase::resetRepositoryMasses() {
  map<string,ParamMap>::const_iterator fit=_parameters.find("mass");
  if(fit==_parameters.end()) 
    throw Exception() << "BLOCK MASS not found in input file"
		      << " can't set masses of SUSY particles"
		      << Exception::runerror;
  ParamMap theMasses = fit->second;
  for(ParamMap::iterator it = theMasses.begin(); it != theMasses.end();){
    long id = it->first;
    double mass = it->second;
    //a negative mass requires an adjustment to the associated mixing matrix 
    if(mass < 0.0) adjustMixingMatrix(id);
    ostringstream cmd;
    tPDPtr part=getParticleData(id);
    if(!part) throw Exception() << "Particle with PDG code " << id 
				<< " not found in SusyBase::"
				<< "resetRepositoryMasses()"
				<< Exception::runerror;
    cmd << "set /Defaults/Particles/" << part->PDGName() 
	<< ":NominalMass " << abs(it->second);
    if(!cmd) 
      throw SetupException() << "SusyBase::resetRepositoryMasses - "
			     << "The stream constructed to reset particle "
			     << "masses is not usuable. " 
			     << Exception::runerror;
    Repository::exec(cmd.str(), cerr);
    //no need to store all masses after this stage
    theMasses.erase(it++);
  }
}

void SusyBase::adjustMixingMatrix(long id) {
  //get correct mixing matrix
  switch(id) {
  case 1000022 :
  case 1000023 :
  case 1000025 :
  case 1000035 : 
  case 1000045 : 
    if(theNMix)
      theNMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The neutralino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  case 1000024 :
  case 1000037 : 
    if(theUMix)
      theUMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The U-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    if(theVMix)
      theVMix->adjustPhase(id);
    else 
      throw SetupException() << "SusyBase::adjustMixingMatrix - "
			     << "The V-Type chargino mixing matrix pointer "
			     << "is null!" << Exception::runerror;
    break;
  default : 
    throw SetupException() << "SusyBase::adjustMixingMatrix - "
			   << "PDG code does not correspond to anything that "
			   << "has a mixing matrix associated with it"
			   << Exception::runerror;
  }
}

void SusyBase::createMixingMatrices() {
  map<string,pair<pair<unsigned int,unsigned int>, vector<MixingElement> > >
    ::const_iterator it;
  for(it=_mixings.begin();it!=_mixings.end();++it) {
    string name=it->first;
    // create the gaugino mixing matrices
    if(name == "nmix" || name == "nmnmix" ){
      createMixingMatrix(theNMix,name,it->second.second,it->second.first);
    }
    else if (name == "umix" ) {
      createMixingMatrix(theUMix,name,it->second.second,it->second.first);
    }
    else if (name == "vmix") {
      createMixingMatrix(theVMix,name,it->second.second,it->second.first);
    }
  }
}

void SusyBase::extractParameters(bool checkmodel) {
  map<string,ParamMap>::const_iterator pit;
  // extract parameters from minpar
  pit=_parameters.find("minpar");
  if(pit==_parameters.end()) 
    throw Exception() << "BLOCK MINPAR not found in " 
		      << "SusyBase::extractParameters()"
		      << Exception::runerror;
  // extract tan beta
  ParamMap::const_iterator it = pit->second.find(3);
  if(it==pit->second.end()) 
    throw Exception() << "Can't find tan beta in BLOCK MINPAR"
		      << Exception::runerror;
  _tanbeta=it->second;
  // extract parameters from hmix
  pit=_parameters.find("hmix");
  if(pit==_parameters.end()) {
    cerr << "BLOCK HMIX not found setting mu to zero\n";
    _mu=0.*GeV;
  }
  else {
    it = pit->second.find(1);
    if(it==pit->second.end()) {
      cerr << "mu not found in BLOCK HMIX setting to zero\n";
    }
    else {
      _mu=it->second*GeV;
    }
  }
  pit = _parameters.find("msoft");
  if( pit == _parameters.end() )
    throw Exception() << "BLOCK MSOFT not found in " 
		      << "SusyBase::extractParameters()"
		      << Exception::runerror;
  it = pit->second.find(1);
  theMone = it->second*GeV;
  it = pit->second.find(2);
  theMtwo = it->second*GeV;
  it = pit->second.find(3);
  theMthree = it->second*GeV;
  if(checkmodel) {
    throw Exception() << "The SusyBase class should not be used as a "
		      << "Model class, use one of the models which inherit"
		      << " from it" << Exception::runerror;
  }
}
