// -*- C++ -*-
//
// RootInterface.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "Herwig++/Interfaces/RootInterface.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

using namespace Herwig;

void RootInterface::init(string filename, string treename){
    f = new TFile(filename.c_str(),"RECREATE"); 
    tree = new TTree(treename.c_str(),treename.c_str());
    chain = new TChain(treename.c_str(),"" );
}

void RootInterface::finish(){
  if (tree->GetFileNumber() !=0) {
    cout << "trying to write and close last file..." << endl;
    gFile->Write();
    gFile->Close();
    cout << "closed the last file" << endl;
    chain->Add("*.root");
    chain->MakeClass();
  } else {
    tree->MakeClass();
    f->Write();
    f->Close();
  }
}

void RootInterface::fill(const string branchname, const double &content){
    map <const string, double> :: const_iterator it = value.find(branchname);

    if ( it == value.end()){
      value[branchname] = content;
      tree->Branch(branchname.c_str(),&value[branchname],(branchname+"/D").c_str() );
    }else{
      value[branchname] = content;
    }
}


void RootInterface::fillarray(const string branchname, const int &i, const double &px,
			      const double &py, const double &pz, const double &E){

  map <const string, TClonesArray*> :: const_iterator it = array.find(branchname);

  if ( it == array.end()){
    array[branchname] = new TClonesArray("TLorentzVector", 500);
    tree->Branch(branchname.c_str(), &array[branchname],32000,0);
  }
  new((*array[branchname])[i]) TLorentzVector(px,py,pz,E);
  
}

void RootInterface::write(){
  /*
    This method is needed, because ::fill cannot assign the values into the tree.
    It would be called several times (as many branches you fill)
  */
  tree->Fill();
  for(map <const string, TClonesArray* >::const_iterator mit = array.begin(); 
      mit!=array.end(); ++mit){
    mit->second->Delete();
  }
}
