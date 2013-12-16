#include "BasicHistogramCollection.h"
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>

using namespace arnold;


inline double deltaphi(double phi1, double phi2){
  double diff=phi1-phi2;
  if(diff<-Constants::pi){
    diff+=(2.0*Constants::pi);
  }
  else if (diff>Constants::pi){
    diff-=(2.0*Constants::pi);
  }
  return diff;
}

double deltaphi(const fastjet::PseudoJet & j1, const fastjet::PseudoJet & j2){
  if (j1.pseudorapidity() > j2.pseudorapidity())
    return deltaphi(j1.phi_std(), j2.phi_std());
  else
    return deltaphi(j2.phi_std(), j1.phi_std());
}

IBPtr BasicHistogramCollection::clone() const {
  return new_ptr(*this);
}

IBPtr BasicHistogramCollection::fullclone() const {
  return new_ptr(*this);
}

void BasicHistogramCollection::fill(const vector<fastjet::PseudoJet> & passed_jets, const vector<fastjet::PseudoJet> & passed_partons, const fastjet::PseudoJet & higgs, const double w){
  if(passed_jets.size()<1 && passed_partons.size()<1){
    cerr<<"\n Warning: BasicHistogramCollection::fill: Size of passed arrays too small! Not filling any histograms for this event.";
    return;
  }


  vector<fastjet::PseudoJet>jets=sorted_by_pt(passed_jets);
  vector<fastjet::PseudoJet>partons=sorted_by_pt(passed_partons);

  size_t i = 0;
  for ( ; i != jets.size(); i++)
    (*histo[96]).addWeighted(i,w);
  (*histo[97]).addWeighted(i-1,w);
  

  if (jets.size() > 0) {
    (*histo[0]).addWeighted(jets[0].perp(),w);
    (*histo[1]).addWeighted(jets[0].eta(),w);
    (*histo[2]).addWeighted(jets[0].rapidity(),w);
    (*histo[3]).addWeighted(jets[0].phi_std(),w);
  }
  if (jets.size() > 1) {
    (*histo[4]).addWeighted(jets[1].perp(),w);
    (*histo[5]).addWeighted(jets[1].eta(),w);
    (*histo[6]).addWeighted(jets[1].rapidity(),w);
    (*histo[7]).addWeighted(jets[1].phi_std(),w);
  }
  if (jets.size() > 2) {
    (*histo[8]).addWeighted(jets[2].perp(),w);
    (*histo[9]).addWeighted(jets[2].eta(),w);
    (*histo[10]).addWeighted(jets[2].rapidity(),w);
    (*histo[11]).addWeighted(jets[2].phi_std(),w);
    (*histo[98]).addWeighted(jets[2].eta(),1.0);
    (*histo[99]).addWeighted(jets[2].eta(),1.0);
  }
  if (jets.size() > 3) {
    (*histo[12]).addWeighted(jets[3].perp(),w);
    (*histo[13]).addWeighted(jets[3].eta(),w);
    (*histo[14]).addWeighted(jets[3].rapidity(),w);
    (*histo[15]).addWeighted(jets[3].phi_std(),w);
  }
  if (jets.size() > 4) {
    (*histo[16]).addWeighted(jets[4].perp(),w);
    (*histo[17]).addWeighted(jets[4].eta(),w);
    (*histo[18]).addWeighted(jets[4].rapidity(),w);
    (*histo[19]).addWeighted(jets[4].phi_std(),w);
  }

  if (jets.size() > 1) {
    (*histo[20]).addWeighted(abs(jets[0].rapidity()-jets[1].rapidity()),w);
    double delEta = abs(jets[0].pseudorapidity()-jets[1].pseudorapidity());
    (*histo[21]).addWeighted(delEta,w);
    double delPhi = deltaphi(jets[0],jets[1]);
    (*histo[22]).addWeighted(delPhi,w);
    double Rjj=sqrt(pow(deltaphi(jets[0],jets[1]),2)+pow(jets[0].pseudorapidity()-jets[1].pseudorapidity(),2));
    (*histo[23]).addWeighted(Rjj,w);
    (*histo[24]).addWeighted((jets[0]+jets[1]).m(),w);

    for (vector<fastjet::PseudoJet>::const_iterator i = jets.begin()+1; i != jets.end(); i++){
      if (i->pseudorapidity()*jets[0].pseudorapidity() > 0){
	(*histo2D[2]).addWeighted(deltaphi(*i,jets[0]),abs(i->pseudorapidity()-jets[0].pseudorapidity()),w);
	break;
      }
    }
    for (vector<fastjet::PseudoJet>::const_iterator i = jets.begin()+1; i != jets.end(); i++){
      if (i->pseudorapidity()*jets[0].pseudorapidity() < 0){
	(*histo2D[3]).addWeighted(deltaphi(*i,jets[0]),abs(i->pseudorapidity()-jets[0].pseudorapidity()),w);
	break;
      }
    }

    (*histo2D[0]).addWeighted(delPhi,delEta,w);
    if ( delEta < 0.4 )
      (*histo[186]).addWeighted(delPhi,w);
    else if ( delEta < 0.8 )
      (*histo[187]).addWeighted(delPhi,w);
    else if ( delEta < 1.2 )
      (*histo[188]).addWeighted(delPhi,w);
    else if ( delEta < 1.6 )
      (*histo[189]).addWeighted(delPhi,w);
    else if ( delEta < 2.0 )
      (*histo[190]).addWeighted(delPhi,w);
    else if ( delEta < 2.4 )
      (*histo[191]).addWeighted(delPhi,w);
    else if ( delEta < 2.8 )
      (*histo[192]).addWeighted(delPhi,w);
    else if ( delEta < 3.2 )
      (*histo[193]).addWeighted(delPhi,w);
    else if ( delEta < 3.6 )
      (*histo[194]).addWeighted(delPhi,w);
    else if ( delEta < 4.0 )
      (*histo[195]).addWeighted(delPhi,w);
    else if ( delEta < 4.4 )
      (*histo[196]).addWeighted(delPhi,w);
    else if ( delEta < 4.8 )
      (*histo[197]).addWeighted(delPhi,w);
    else if ( delEta < 5.2 )
      (*histo[198]).addWeighted(delPhi,w);
    else if ( delEta < 5.6 )
      (*histo[199]).addWeighted(delPhi,w);
    else if ( delEta < 6.0 )
      (*histo[200]).addWeighted(delPhi,w);
    else if ( delEta < 6.4 )
      (*histo[201]).addWeighted(delPhi,w);
    else if ( delEta < 6.8 )
      (*histo[202]).addWeighted(delPhi,w);
    else if ( delEta < 7.2 )
      (*histo[203]).addWeighted(delPhi,w);
    else if ( delEta < 7.6 )
      (*histo[204]).addWeighted(delPhi,w);
    else if ( delEta < 8.0 )
      (*histo[205]).addWeighted(delPhi,w);
    else if ( delEta < 8.4 )
      (*histo[206]).addWeighted(delPhi,w);
    else if ( delEta < 8.8 )
      (*histo[207]).addWeighted(delPhi,w);
    else if ( delEta < 9.2 )
      (*histo[208]).addWeighted(delPhi,w);
    else if ( delEta < 9.6 )
      (*histo[209]).addWeighted(delPhi,w);
    else if ( delEta < 10.0 )
      (*histo[210]).addWeighted(delPhi,w);
  }

  if (jets.size() > 2) {
    (*histo[25]).addWeighted(abs(jets[0].rapidity()-jets[2].rapidity()),w);
    double delEta = abs(jets[0].pseudorapidity()-jets[2].pseudorapidity());
    (*histo[26]).addWeighted(delEta,w);
    double delPhi = deltaphi(jets[0],jets[2]);
    (*histo[27]).addWeighted(delPhi,w);
    double Rjj=sqrt(pow(deltaphi(jets[0],jets[2]),2)+pow(jets[0].pseudorapidity()-jets[2].pseudorapidity(),2));
    (*histo[28]).addWeighted(Rjj,w);
    (*histo[29]).addWeighted((jets[0]+jets[2]).m(),w);


    for (vector<fastjet::PseudoJet>::const_iterator i = jets.begin(); i != jets.end(); i++){
      if (*i == jets[1]) continue;
      if (i->pseudorapidity()*jets[1].pseudorapidity() > 0){
	(*histo2D[4]).addWeighted(deltaphi(*i,jets[1]),abs(i->pseudorapidity()-jets[1].pseudorapidity()),w);
	break;
      }
    }
    for (vector<fastjet::PseudoJet>::const_iterator i = jets.begin()+1; i != jets.end(); i++){
      if (*i == jets[1]) continue;
      if (i->pseudorapidity()*jets[1].pseudorapidity() < 0){
	(*histo2D[5]).addWeighted(deltaphi(*i,jets[1]),abs(i->pseudorapidity()-jets[1].pseudorapidity()),w);
	break;
      }
    }
    
    double delPhi12 = deltaphi(jets[0],jets[1]);

    (*histo2D[6]).addWeighted(delPhi12,jets[2].perp(),w);
    (*histo2D[7]).addWeighted(delPhi,jets[2].perp(),w);

    (*histo2D[1]).addWeighted(delPhi,delEta,w);
    if ( delEta < 0.4 )
      (*histo[211]).addWeighted(delPhi,w);
    else if ( delEta < 0.8 )
      (*histo[212]).addWeighted(delPhi,w);
    else if ( delEta < 1.2 )
      (*histo[213]).addWeighted(delPhi,w);
    else if ( delEta < 1.6 )
      (*histo[214]).addWeighted(delPhi,w);
    else if ( delEta < 2.0 )
      (*histo[215]).addWeighted(delPhi,w);
    else if ( delEta < 2.4 )
      (*histo[216]).addWeighted(delPhi,w);
    else if ( delEta < 2.8 )
      (*histo[217]).addWeighted(delPhi,w);
    else if ( delEta < 3.2 )
      (*histo[218]).addWeighted(delPhi,w);
    else if ( delEta < 3.6 )
      (*histo[219]).addWeighted(delPhi,w);
    else if ( delEta < 4.0 )
      (*histo[220]).addWeighted(delPhi,w);
    else if ( delEta < 4.4 )
      (*histo[221]).addWeighted(delPhi,w);
    else if ( delEta < 4.8 )
      (*histo[222]).addWeighted(delPhi,w);
    else if ( delEta < 5.2 )
      (*histo[223]).addWeighted(delPhi,w);
    else if ( delEta < 5.6 )
      (*histo[224]).addWeighted(delPhi,w);
    else if ( delEta < 6.0 )
      (*histo[225]).addWeighted(delPhi,w);
    else if ( delEta < 6.4 )
      (*histo[226]).addWeighted(delPhi,w);
    else if ( delEta < 6.8 )
      (*histo[227]).addWeighted(delPhi,w);
    else if ( delEta < 7.2 )
      (*histo[228]).addWeighted(delPhi,w);
    else if ( delEta < 7.6 )
      (*histo[229]).addWeighted(delPhi,w);
    else if ( delEta < 8.0 )
      (*histo[230]).addWeighted(delPhi,w);
    else if ( delEta < 8.4 )
      (*histo[231]).addWeighted(delPhi,w);
    else if ( delEta < 8.8 )
      (*histo[232]).addWeighted(delPhi,w);
    else if ( delEta < 9.2 )
      (*histo[233]).addWeighted(delPhi,w);
    else if ( delEta < 9.6 )
      (*histo[234]).addWeighted(delPhi,w);
    else if ( delEta < 10.0 )
      (*histo[235]).addWeighted(delPhi,w);
  }

  if (jets.size() > 3) {
    (*histo[30]).addWeighted(abs(jets[0].rapidity()-jets[3].rapidity()),w);
    (*histo[31]).addWeighted(abs(jets[0].pseudorapidity()-jets[3].pseudorapidity()),w);
    (*histo[32]).addWeighted(deltaphi(jets[0],jets[3]),w);
    double Rjj=sqrt(pow(deltaphi(jets[0],jets[3]),2)+pow(jets[0].pseudorapidity()-jets[3].pseudorapidity(),2));
    (*histo[33]).addWeighted(Rjj,w);
    (*histo[34]).addWeighted((jets[0]+jets[3]).m(),w);
  }

  if (jets.size() > 2) {
    (*histo[35]).addWeighted(abs(jets[1].rapidity()-jets[2].rapidity()),w);
    (*histo[36]).addWeighted(abs(jets[1].pseudorapidity()-jets[2].pseudorapidity()),w);
    (*histo[37]).addWeighted(deltaphi(jets[1],jets[2]),w);
    double Rjj=sqrt(pow(deltaphi(jets[1],jets[2]),2)+pow(jets[1].pseudorapidity()-jets[2].pseudorapidity(),2));
    (*histo[38]).addWeighted(Rjj,w);
    (*histo[39]).addWeighted((jets[1]+jets[2]).m(),w);
  }

  if (jets.size() > 3) {
    (*histo[40]).addWeighted(abs(jets[1].rapidity()-jets[3].rapidity()),w);
    (*histo[41]).addWeighted(abs(jets[1].pseudorapidity()-jets[3].pseudorapidity()),w);
    (*histo[42]).addWeighted(deltaphi(jets[1],jets[3]),w);
    double Rjj=sqrt(pow(deltaphi(jets[1],jets[3]),2)+pow(jets[1].pseudorapidity()-jets[3].pseudorapidity(),2));
    (*histo[43]).addWeighted(Rjj,w);
    (*histo[44]).addWeighted((jets[1]+jets[3]).m(),w);
  }

  if (jets.size() > 3) {
    (*histo[45]).addWeighted(abs(jets[2].rapidity()-jets[3].rapidity()),w);
    (*histo[46]).addWeighted(abs(jets[2].pseudorapidity()-jets[3].pseudorapidity()),w);
    (*histo[47]).addWeighted(deltaphi(jets[2],jets[3]),w);
    double Rjj=sqrt(pow(deltaphi(jets[2],jets[3]),2)+pow(jets[2].pseudorapidity()-jets[3].pseudorapidity(),2));
    (*histo[48]).addWeighted(Rjj,w);
    (*histo[49]).addWeighted((jets[2]+jets[3]).m(),w);
  }

  if (jets.size() > 0) {
    (*histo[50]).addWeighted(abs(higgs.rapidity()-jets[0].rapidity()),w);
    (*histo[51]).addWeighted(deltaphi(higgs,jets[0]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,jets[0]),2)+pow(higgs.pseudorapidity()-jets[0].pseudorapidity(),2));
    (*histo[52]).addWeighted(Rjj,w);
    (*histo[53]).addWeighted((higgs+jets[0]).m(),w);
  }

  if (jets.size() > 1) {
    (*histo[54]).addWeighted(abs(higgs.rapidity()-jets[1].rapidity()),w);
    (*histo[55]).addWeighted(deltaphi(higgs,jets[1]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,jets[1]),2)+pow(higgs.pseudorapidity()-jets[1].pseudorapidity(),2));
    (*histo[56]).addWeighted(Rjj,w);
    (*histo[57]).addWeighted((higgs+jets[1]).m(),w);
  }

  if (jets.size() > 2) {
    (*histo[58]).addWeighted(abs(higgs.rapidity()-jets[2].rapidity()),w);
    (*histo[59]).addWeighted(deltaphi(higgs,jets[2]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,jets[2]),2)+pow(higgs.pseudorapidity()-jets[2].pseudorapidity(),2));
    (*histo[60]).addWeighted(Rjj,w);
    (*histo[61]).addWeighted((higgs+jets[2]).m(),w);
  }

  if (jets.size() > 3) {
    (*histo[62]).addWeighted(abs(higgs.rapidity()-jets[3].rapidity()),w);
    (*histo[63]).addWeighted(deltaphi(higgs,jets[3]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,jets[3]),2)+pow(higgs.pseudorapidity()-jets[3].pseudorapidity(),2));
    (*histo[64]).addWeighted(Rjj,w);
    (*histo[65]).addWeighted((higgs+jets[3]).m(),w);
  }

  if (jets.size() > 2) {
    (*histo[66]).addWeighted(jets[2].perp()-(jets[0].perp()+jets[1].perp())/2.,w);
    (*histo[67]).addWeighted(jets[2].rapidity()-(jets[0].rapidity()+jets[1].rapidity())/2.,w);
    (*histo[68]).addWeighted(jets[2].pseudorapidity()-(jets[0].pseudorapidity()+jets[1].pseudorapidity())/2.,w);
    (*histo[69]).addWeighted(abs(jets[2].phi_std() - (abs(jets[0].phi_std())+abs(jets[1].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = jets[0] + jets[1];
    double Rjj=sqrt(pow(deltaphi(jets[2],combined),2)+pow(jets[2].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[70]).addWeighted(Rjj,w);
  }

  if (jets.size() > 3) {
    (*histo[71]).addWeighted(jets[3].perp()-(jets[0].perp()+jets[1].perp())/2.,w);
    (*histo[72]).addWeighted(jets[3].rapidity()-(jets[0].rapidity()+jets[1].rapidity())/2.,w);
    (*histo[73]).addWeighted(jets[3].pseudorapidity()-(jets[0].pseudorapidity()+jets[1].pseudorapidity())/2.,w);
    (*histo[74]).addWeighted(abs(jets[3].phi_std() - (abs(jets[0].phi_std())+abs(jets[1].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = jets[0] + jets[1];
    double Rjj=sqrt(pow(deltaphi(jets[3],combined),2)+pow(jets[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[75]).addWeighted(Rjj,w);
  }

  if (jets.size() > 3) {
    (*histo[76]).addWeighted(jets[3].perp()-(jets[0].perp()+jets[2].perp())/2.,w);
    (*histo[77]).addWeighted(jets[3].rapidity()-(jets[0].rapidity()+jets[2].rapidity())/2.,w);
    (*histo[78]).addWeighted(jets[3].pseudorapidity()-(jets[0].pseudorapidity()+jets[2].pseudorapidity())/2.,w);
    (*histo[79]).addWeighted(abs(jets[3].phi_std() - (abs(jets[0].phi_std())+abs(jets[2].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = jets[0] + jets[2];
    double Rjj=sqrt(pow(deltaphi(jets[3],combined),2)+pow(jets[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[80]).addWeighted(Rjj,w);
  }

  if (jets.size() > 3) {
    (*histo[81]).addWeighted(jets[3].perp()-(jets[1].perp()+jets[2].perp())/2.,w);
    (*histo[82]).addWeighted(jets[3].rapidity()-(jets[1].rapidity()+jets[2].rapidity())/2.,w);
    (*histo[83]).addWeighted(jets[3].pseudorapidity()-(jets[1].pseudorapidity()+jets[2].pseudorapidity())/2.,w);
    (*histo[84]).addWeighted(abs(jets[3].phi_std() - (abs(jets[1].phi_std())+abs(jets[2].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = jets[1] + jets[2];
    double Rjj=sqrt(pow(deltaphi(jets[3],combined),2)+pow(jets[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[85]).addWeighted(Rjj,w);
  }

  if (jets.size() > 2) {
    double ystar = jets[2].rapidity() - 0.5 * (jets[0].rapidity()+jets[1].rapidity());
    (*histo[86]).addWeighted(ystar,w);
    if ( abs(jets[0].rapidity()-jets[1].rapidity()) > 0.000001 )
      (*histo[87]).addWeighted(ystar/(abs(jets[0].rapidity()-jets[1].rapidity())),w);
  }
  
  unsigned int n1min=-1;
  unsigned int n2min=-1;
  unsigned int n3min=-1;
  unsigned int n4min=-1;
  double  Rp1jn=100;
  double  Rp2jn=100;
  double  Rp3jn=100;
  double  Rp4jn=100;
  
  if (partons.size() > 0){
    for(unsigned int i = 0; i < jets.size(); i++){
      if ( sqrt(pow(deltaphi(partons[0],jets[i]),2)+pow(partons[0].rapidity()-jets[i].rapidity(),2)) < Rp1jn ){
	Rp1jn = sqrt(pow(deltaphi(partons[0],jets[i]),2)+pow(partons[0].rapidity()-jets[i].rapidity(),2));
	n1min = i;
      }
    }
    if ( Rp1jn > 1.0 )
      n1min = -1;
    (*histo[88]).addWeighted(n1min,w);
  }  

  if (partons.size() > 1){
    for(unsigned int i = 0; i < jets.size(); i++){
      if ( sqrt(pow(deltaphi(partons[1],jets[i]),2)+pow(partons[1].rapidity()-jets[i].rapidity(),2)) < Rp1jn ){
	Rp2jn = sqrt(pow(deltaphi(partons[1],jets[i]),2)+pow(partons[1].rapidity()-jets[i].rapidity(),2));
	n2min = i;
      }
    }
    if ( Rp2jn > 1.0 )
      n2min = -1;
    (*histo[89]).addWeighted(n2min,w);
  }  

  if (partons.size() > 2){
    for(unsigned int i = 0; i < jets.size(); i++){
      if ( sqrt(pow(deltaphi(partons[2],jets[i]),2)+pow(partons[2].rapidity()-jets[i].rapidity(),2)) < Rp1jn ){
	Rp3jn = sqrt(pow(deltaphi(partons[2],jets[i]),2)+pow(partons[2].rapidity()-jets[i].rapidity(),2));
	n3min = i;
      }
    }
    if ( Rp3jn > 1.0 )
      n3min = -1;
    (*histo[90]).addWeighted(n3min,w);
  }  

  if (partons.size() > 3){
    for(unsigned int i = 0; i < jets.size(); i++){
      if ( sqrt(pow(deltaphi(partons[3],jets[i]),2)+pow(partons[3].rapidity()-jets[i].rapidity(),2)) < Rp1jn ){
	Rp4jn = sqrt(pow(deltaphi(partons[3],jets[i]),2)+pow(partons[3].rapidity()-jets[i].rapidity(),2));
	n4min = i;
      }
    }
    if ( Rp4jn > 1.0 )
      n4min = -1;
    (*histo[91]).addWeighted(n4min,w);
  }

  (*histo[92]).addWeighted(higgs.perp(),w);
  (*histo[93]).addWeighted(higgs.eta(),w);
  (*histo[94]).addWeighted(higgs.rapidity(),w);
  (*histo[95]).addWeighted(higgs.phi_std(),w);
  
  if (partons.size() > 0) {
    (*histo[100]).addWeighted(partons[0].perp(),w);
    (*histo[101]).addWeighted(partons[0].eta(),w);
    (*histo[102]).addWeighted(partons[0].rapidity(),w);
    (*histo[103]).addWeighted(partons[0].phi_std(),w);
  }
  if (partons.size() > 1) {
    (*histo[104]).addWeighted(partons[1].perp(),w);
    (*histo[105]).addWeighted(partons[1].eta(),w);
    (*histo[106]).addWeighted(partons[1].rapidity(),w);
    (*histo[107]).addWeighted(partons[1].phi_std(),w);
  }
  if (partons.size() > 2) {
    (*histo[108]).addWeighted(partons[2].perp(),w);
    (*histo[109]).addWeighted(partons[2].eta(),w);
    (*histo[110]).addWeighted(partons[2].rapidity(),w);
    (*histo[111]).addWeighted(partons[2].phi_std(),w);
  }
  if (partons.size() > 3) {
    (*histo[112]).addWeighted(partons[3].perp(),w);
    (*histo[113]).addWeighted(partons[3].eta(),w);
    (*histo[114]).addWeighted(partons[3].rapidity(),w);
    (*histo[115]).addWeighted(partons[3].phi_std(),w);
  }
  if (partons.size() > 4) {
    (*histo[116]).addWeighted(partons[4].perp(),w);
    (*histo[117]).addWeighted(partons[4].eta(),w);
    (*histo[118]).addWeighted(partons[4].rapidity(),w);
    (*histo[119]).addWeighted(partons[4].phi_std(),w);
  }

  if (partons.size() > 1) {
    (*histo[120]).addWeighted(abs(partons[0].rapidity()-partons[1].rapidity()),w);
    (*histo[121]).addWeighted(abs(partons[0].pseudorapidity()-partons[1].pseudorapidity()),w);
    (*histo[122]).addWeighted(deltaphi(partons[0],partons[1]),w);
    double Rjj=sqrt(pow(deltaphi(partons[0],partons[1]),2)+pow(partons[0].pseudorapidity()-partons[1].pseudorapidity(),2));
    (*histo[123]).addWeighted(Rjj,w);
    (*histo[124]).addWeighted((partons[0]+partons[1]).m(),w);
  }

  if (partons.size() > 2) {
    (*histo[125]).addWeighted(abs(partons[0].rapidity()-partons[2].rapidity()),w);
    (*histo[126]).addWeighted(abs(partons[0].pseudorapidity()-partons[2].pseudorapidity()),w);
    (*histo[127]).addWeighted(deltaphi(partons[0],partons[2]),w);
    double Rjj=sqrt(pow(deltaphi(partons[0],partons[2]),2)+pow(partons[0].pseudorapidity()-partons[2].pseudorapidity(),2));
    (*histo[128]).addWeighted(Rjj,w);
    (*histo[129]).addWeighted((partons[0]+partons[2]).m(),w);
  }

  if (partons.size() > 3) {
    (*histo[130]).addWeighted(abs(partons[0].rapidity()-partons[3].rapidity()),w);
    (*histo[131]).addWeighted(abs(partons[0].pseudorapidity()-partons[3].pseudorapidity()),w);
    (*histo[132]).addWeighted(deltaphi(partons[0],partons[3]),w);
    double Rjj=sqrt(pow(deltaphi(partons[0],partons[3]),2)+pow(partons[0].pseudorapidity()-partons[3].pseudorapidity(),2));
    (*histo[133]).addWeighted(Rjj,w);
    (*histo[134]).addWeighted((partons[0]+partons[3]).m(),w);
  }

  if (partons.size() > 2) {
    (*histo[135]).addWeighted(abs(partons[1].rapidity()-partons[2].rapidity()),w);
    (*histo[136]).addWeighted(abs(partons[1].pseudorapidity()-partons[2].pseudorapidity()),w);
    (*histo[137]).addWeighted(deltaphi(partons[1],partons[2]),w);
    double Rjj=sqrt(pow(deltaphi(partons[1],partons[2]),2)+pow(partons[1].pseudorapidity()-partons[2].pseudorapidity(),2));
    (*histo[138]).addWeighted(Rjj,w);
    (*histo[139]).addWeighted((partons[1]+partons[2]).m(),w);
  }

  if (partons.size() > 3) {
    (*histo[140]).addWeighted(abs(partons[1].rapidity()-partons[3].rapidity()),w);
    (*histo[141]).addWeighted(abs(partons[1].pseudorapidity()-partons[3].pseudorapidity()),w);
    (*histo[142]).addWeighted(deltaphi(partons[1],partons[3]),w);
    double Rjj=sqrt(pow(deltaphi(partons[1],partons[3]),2)+pow(partons[1].pseudorapidity()-partons[3].pseudorapidity(),2));
    (*histo[143]).addWeighted(Rjj,w);
    (*histo[144]).addWeighted((partons[1]+partons[3]).m(),w);
  }

  if (partons.size() > 3) {
    (*histo[145]).addWeighted(abs(partons[2].rapidity()-partons[3].rapidity()),w);
    (*histo[146]).addWeighted(abs(partons[2].pseudorapidity()-partons[3].pseudorapidity()),w);
    (*histo[147]).addWeighted(deltaphi(partons[2],partons[3]),w);
    double Rjj=sqrt(pow(deltaphi(partons[2],partons[3]),2)+pow(partons[2].pseudorapidity()-partons[3].pseudorapidity(),2));
    (*histo[148]).addWeighted(Rjj,w);
    (*histo[149]).addWeighted((partons[2]+partons[3]).m(),w);
  }

  if (partons.size() > 0) {
    (*histo[150]).addWeighted(abs(higgs.rapidity()-partons[0].rapidity()),w);
    (*histo[151]).addWeighted(deltaphi(higgs,partons[0]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,partons[0]),2)+pow(higgs.pseudorapidity()-partons[0].pseudorapidity(),2));
    (*histo[152]).addWeighted(Rjj,w);
    (*histo[153]).addWeighted((higgs+partons[0]).m(),w);
  }

  if (partons.size() > 1) {
    (*histo[154]).addWeighted(abs(higgs.rapidity()-partons[1].rapidity()),w);
    (*histo[155]).addWeighted(deltaphi(higgs,partons[1]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,partons[1]),2)+pow(higgs.pseudorapidity()-partons[1].pseudorapidity(),2));
    (*histo[156]).addWeighted(Rjj,w);
    (*histo[157]).addWeighted((higgs+partons[1]).m(),w);
  }

  if (partons.size() > 2) {
    (*histo[158]).addWeighted(abs(higgs.rapidity()-partons[2].rapidity()),w);
    (*histo[159]).addWeighted(deltaphi(higgs,partons[2]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,partons[2]),2)+pow(higgs.pseudorapidity()-partons[2].pseudorapidity(),2));
    (*histo[160]).addWeighted(Rjj,w);
    (*histo[161]).addWeighted((higgs+partons[2]).m(),w);
  }

  if (partons.size() > 3) {
    (*histo[162]).addWeighted(abs(higgs.rapidity()-partons[3].rapidity()),w);
    (*histo[163]).addWeighted(deltaphi(higgs,partons[3]),w);
    double Rjj=sqrt(pow(deltaphi(higgs,partons[3]),2)+pow(higgs.pseudorapidity()-partons[3].pseudorapidity(),2));
    (*histo[164]).addWeighted(Rjj,w);
    (*histo[165]).addWeighted((higgs+partons[3]).m(),w);
  }

  if (partons.size() > 2) {
    (*histo[166]).addWeighted(partons[2].perp()-(partons[0].perp()+partons[1].perp())/2.,w);
    (*histo[167]).addWeighted(partons[2].rapidity()-(partons[0].rapidity()+partons[1].rapidity())/2.,w);
    (*histo[168]).addWeighted(partons[2].pseudorapidity()-(partons[0].pseudorapidity()+partons[1].pseudorapidity())/2.,w);
    (*histo[169]).addWeighted(abs(partons[2].phi_std() - (abs(partons[0].phi_std())+abs(partons[1].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = partons[0] + partons[1];
    double Rjj=sqrt(pow(deltaphi(partons[2],combined),2)+pow(partons[2].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[170]).addWeighted(Rjj,w);
  }

  if (partons.size() > 3) {
    (*histo[171]).addWeighted(partons[3].perp()-(partons[0].perp()+partons[1].perp())/2.,w);
    (*histo[172]).addWeighted(partons[3].rapidity()-(partons[0].rapidity()+partons[1].rapidity())/2.,w);
    (*histo[173]).addWeighted(partons[3].pseudorapidity()-(partons[0].pseudorapidity()+partons[1].pseudorapidity())/2.,w);
    (*histo[174]).addWeighted(abs(partons[3].phi_std() - (abs(partons[0].phi_std())+abs(partons[1].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = partons[0] + partons[1];
    double Rjj=sqrt(pow(deltaphi(partons[3],combined),2)+pow(partons[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[175]).addWeighted(Rjj,w);
  }

  if (partons.size() > 3) {
    (*histo[176]).addWeighted(partons[3].perp()-(partons[0].perp()+partons[2].perp())/2.,w);
    (*histo[177]).addWeighted(partons[3].rapidity()-(partons[0].rapidity()+partons[2].rapidity())/2.,w);
    (*histo[178]).addWeighted(partons[3].pseudorapidity()-(partons[0].pseudorapidity()+partons[2].pseudorapidity())/2.,w);
    (*histo[179]).addWeighted(abs(partons[3].phi_std() - (abs(partons[0].phi_std())+abs(partons[2].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = partons[0] + partons[2];
    double Rjj=sqrt(pow(deltaphi(partons[3],combined),2)+pow(partons[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[180]).addWeighted(Rjj,w);
  }

  if (partons.size() > 3) {
    (*histo[181]).addWeighted(partons[3].perp()-(partons[1].perp()+partons[2].perp())/2.,w);
    (*histo[182]).addWeighted(partons[3].rapidity()-(partons[1].rapidity()+partons[2].rapidity())/2.,w);
    (*histo[183]).addWeighted(partons[3].pseudorapidity()-(partons[1].pseudorapidity()+partons[2].pseudorapidity())/2.,w);
    (*histo[184]).addWeighted(abs(partons[3].phi_std() - (abs(partons[1].phi_std())+abs(partons[2].phi_std()))/2 ),w);
    fastjet::PseudoJet combined = partons[1] + partons[2];
    double Rjj=sqrt(pow(deltaphi(partons[3],combined),2)+pow(partons[3].pseudorapidity()-combined.pseudorapidity(),2));
    (*histo[185]).addWeighted(Rjj,w);
  }

  


  // if (jets.size() > 0){
  //   (*histo[1]).addWeighted(jets[0].perp(),w);
  //   (*histo[6]).addWeighted(jets[0].eta(),w);
  //   (*histo[28]).addWeighted(jets[0].rapidity(),w);
  // }

  // if (partons.size() > 0) {
  //   (*histo[11]).addWeighted(partons[0].perp(),w);
  //   (*histo[16]).addWeighted(partons[0].rapidity(),w);
  // }

  // if (jets.size() > 1) {
  //   double Rjj=sqrt(pow(deltaphi(jets[0].phi_std(),jets[1].phi_std()),2)+pow(jets[0].rapidity()-jets[1].rapidity(),2));
  //   (*histo[0]).addWeighted(Rjj,w);
  //   (*histo[2]).addWeighted(jets[1].perp(),w);
  //   (*histo[3]).addWeighted(abs(jets[0].eta()-jets[1].eta()),w);
  //   (*histo[30]).addWeighted(abs(jets[0].rapidity()-jets[1].rapidity()),w);    
  //   fastjet::PseudoJet sumj = jets[0];
  //   sumj+=jets[1];
  //   (*histo[4]).addWeighted(sumj.m(),w);
  //   if(jets[0].rapidity()>jets[1].rapidity()){
  //     (*histo[5]).addWeighted(deltaphi(jets[1].phi_std(),jets[0].phi_std()),w);
  //   }
  //   else{
  //     (*histo[5]).addWeighted(deltaphi(jets[0].phi_std(),jets[1].phi_std()),w);
  //   } 
  //   if(abs(deltaphi(jets[1].phi_std(),jets[0].phi_std()))<0.25*Constants::pi){Nj_delphi_0to25+=1;}
  //   else if(abs(deltaphi(jets[1].phi_std(),jets[0].phi_std()))<0.75*Constants::pi){Nj_delphi_25to75+=1;}
  //   else if(abs(deltaphi(jets[1].phi_std(),jets[0].phi_std()))<Constants::pi){Nj_delphi_75to100+=1;}	
  //   A_phi_jets=(double)(Nj_delphi_0to25 -Nj_delphi_25to75 + Nj_delphi_75to100)/(Nj_delphi_0to25 + Nj_delphi_25to75 + Nj_delphi_75to100);
  //   (*histo[7]).addWeighted(jets[1].eta(),w);
  //   (*histo[29]).addWeighted(jets[1].rapidity(),w);
  //   (*histo[26]).addWeighted(jets[0].perp()/jets[1].perp(),w);

  // }
  // if (partons.size() > 1) {
  //   double Rpp=sqrt(pow(deltaphi(partons[0].phi_std(),partons[1].phi_std()),2)+pow(partons[0].rapidity()-partons[1].rapidity(),2));  
  //   if(abs(deltaphi(partons[1].phi_std(),partons[0].phi_std()))<0.25*Constants::pi){Np_delphi_0to25+=1;}
  //   else if(abs(deltaphi(partons[1].phi_std(),partons[0].phi_std()))<0.75*Constants::pi){Np_delphi_25to75+=1;}
  //   else if(abs(deltaphi(partons[1].phi_std(),partons[0].phi_std()))<Constants::pi){Np_delphi_75to100+=1;}	
    
  //   A_phi_partons=(double)(Np_delphi_0to25 -Np_delphi_25to75 + Np_delphi_75to100)/(Np_delphi_0to25 + Np_delphi_25to75 + Np_delphi_75to100);
    
  //   (*histo[10]).addWeighted(Rpp,w);
  //   (*histo[12]).addWeighted(partons[1].perp(),w);
  //   (*histo[13]).addWeighted(abs(partons[0].rapidity()-partons[1].rapidity()),w);
    
  //   fastjet::PseudoJet sump = partons[0];
  //   sump+=partons[1];
  //   (*histo[14]).addWeighted(sump.m(),w);
  //   if(partons[0].rapidity()>partons[1].rapidity()){
  //     (*histo[15]).addWeighted(deltaphi(partons[1].phi_std(),partons[0].phi_std()),w);
  //   } 
  //   else{
  //     (*histo[15]).addWeighted(deltaphi(partons[0].phi_std(),partons[1].phi_std()),w);
  //   }
  //   (*histo[17]).addWeighted(partons[1].rapidity(),w);
  // }
  // if (jets.size() > 1 && partons.size() > 1) {
  //   double R11=sqrt(pow(deltaphi(jets[0].phi_std(),partons[0].phi_std()),2)+pow(jets[0].rapidity()-partons[0].rapidity(),2));
  //   double R22=sqrt(pow(deltaphi(jets[1].phi_std(),partons[1].phi_std()),2)+pow(jets[1].rapidity()-partons[1].rapidity(),2));
  //   double R21=sqrt(pow(deltaphi(jets[1].phi_std(),partons[0].phi_std()),2)+pow(jets[1].rapidity()-partons[0].rapidity(),2));
  //   double R12=sqrt(pow(deltaphi(jets[0].phi_std(),partons[1].phi_std()),2)+pow(jets[0].rapidity()-partons[1].rapidity(),2));
  //   (*histo[8]).addWeighted(min(R11,R12),w);
  //   (*histo[9]).addWeighted(min(R21,R22),w);
    
  
  //   (*histo[18]).addWeighted(jets[0].perp()-partons[0].perp(),w);
  //   (*histo[19]).addWeighted(jets[1].perp()-partons[1].perp(),w);
  //   (*histo[20]).addWeighted(jets[0].rapidity()-partons[0].rapidity(),w);
  //   (*histo[21]).addWeighted(jets[1].rapidity()-partons[1].rapidity(),w);
  //   (*histo[22]).addWeighted(deltaphi(jets[0].phi_std(),partons[0].phi_std()),w);
  //   (*histo[23]).addWeighted(deltaphi(jets[1].phi_std(),partons[1].phi_std()),w);
  //   (*histo[27]).addWeighted(partons[0].perp()/partons[1].perp(),w);

  // }  
  

  // if (jets.size() > 2) {
  //   (*histo[31]).addWeighted(jets[2].perp(),w);
  //   (*histo[33]).addWeighted(jets[2].eta(),w);
  //   (*histo[35]).addWeighted(jets[2].rapidity(),w);
  //   double ystar = jets[2].rapidity() - 0.5 * (jets[0].rapidity()+jets[1].rapidity());
  //   (*histo[37]).addWeighted(ystar,w);
  //   if ( abs(jets[0].rapidity()-jets[1].rapidity()) > 0.000001 )
  //     (*histo[38]).addWeighted(ystar/(abs(jets[0].rapidity()-jets[1].rapidity())),w);
  // }
  // if (jets.size() > 3) {
  //   (*histo[32]).addWeighted(jets[3].perp(),w);
  //   (*histo[34]).addWeighted(jets[3].eta(),w);
  //   (*histo[36]).addWeighted(jets[3].rapidity(),w);
  // }

  // unsigned int n=-1;
  // unsigned int m=-1;
  // double  Rp1jn=100;
  // double  Rp2jm=100;
  

  // do{
  //   if(n==jets.size()){
  //     n=-1;
  //     break;
  //   }
  //   n++;
  //   Rp1jn=sqrt(pow(deltaphi(partons[0].phi_std(),jets[n].phi_std()),2)+pow(partons[0].rapidity()-jets[n].rapidity(),2));
    
  // }while(Rp1jn>1.0);
  // do{
  //   if(m==jets.size()){
  //     m=-1;
  //     break;
  //   }
  //   m++;
  //   Rp2jm=sqrt(pow(deltaphi(partons[1].phi_std(),jets[m].phi_std()),2)+pow(partons[1].rapidity()-jets[m].rapidity(),2));
  // }while(Rp2jm>1.0);
  

  // (*histo[24]).addWeighted(n+1,w);
  // (*histo[25]).addWeighted(m+1,w);

  // if(n==0 && m==1){          //if p1=j1 & p2=j2
  //   correlarray[0][0]+=1;
  // } 
  // else if(n==1 && m==0){     //if p1=j2 & p2!=j1
  //   correlarray[0][2]+=1;
  // } 
  // else if(n==0 && m!=1){     //if p1=j1 & p2!=j2
  //   correlarray[0][1]+=1;
  // } 
  // else if(n!=0 && m==1){     //if p1!=j1 & p2=j2
  //   correlarray[1][0]+=1;
  // } 
  // else if((n!=0 && m!=1)&&(n!=1 && m!=0)){ //if p1/2!=j1 & p1/2!=j2
  //   correlarray[1][1]+=1;
  // }
  

}



BasicHistogramCollection::BasicHistogramCollection()
  : theNBins(100), doSmearing(true){

  int n=100;

  histo[0]= arHistogramPtr::Create(arHistogram(0,500,n,doSmearing));   //"pt_j1.dat"
  histo[1]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j1.dat"
  histo[2]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j1.dat"
  histo[3]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j1.dat"

  histo[4]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_j2.dat"
  histo[5]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j2.dat"
  histo[6]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j2.dat"
  histo[7]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j2.dat"

  histo[8]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_j3.dat"
  histo[9]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j3.dat"
  histo[10]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j3.dat"
  histo[11]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j3.dat"

  histo[12]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_j4.dat"
  histo[13]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j4.dat"
  histo[14]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j4.dat"
  histo[15]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j4.dat"

  histo[16]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_j5.dat"
  histo[17]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j5.dat"
  histo[18]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j5.dat"
  histo[19]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j5.dat"

  histo[20]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j1j2.dat"
  histo[21]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j1j2.dat"
  histo[22]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j1j2.dat"
  histo[23]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j1j2.dat"
  histo[24]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j1j2.dat"

  histo[25]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j1j3.dat"
  histo[26]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j1j3.dat"
  histo[27]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j1j3.dat"
  histo[28]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j1j3.dat"
  histo[29]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j1j3.dat"

  histo[30]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j1j4.dat"
  histo[31]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j1j4.dat"
  histo[32]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j1j4.dat"
  histo[33]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j1j4.dat"
  histo[34]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j1j4.dat"

  histo[35]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j2j3.dat"
  histo[36]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j2j3.dat"
  histo[37]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j2j3.dat"
  histo[38]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j2j3.dat"
  histo[39]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j2j3.dat"

  histo[40]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j2j4.dat"
  histo[41]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j2j4.dat"
  histo[42]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j2j4.dat"
  histo[43]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j2j4.dat"
  histo[44]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j2j3.dat"

  histo[45]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_j3j4.dat"
  histo[46]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_j3j4.dat"
  histo[47]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_j3j4.dat"
  histo[48]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_j3j4.dat"
  histo[49]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_j3j4.dat"

  histo[50]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hj1.dat"
  histo[51]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hj1.dat"
  histo[52]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hj1.dat"
  histo[53]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hj1.dat"

  histo[54]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hj2.dat"
  histo[55]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hj2.dat"
  histo[56]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hj2.dat"
  histo[57]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hj2.dat"

  histo[58]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hj3.dat"
  histo[59]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hj3.dat"
  histo[60]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hj3.dat"
  histo[61]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hj3.dat"

  histo[62]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hj4.dat"
  histo[63]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hj4.dat"
  histo[64]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hj4.dat"
  histo[65]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hj4.dat"

  histo[66]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_j3VSj1j2.dat"
  histo[67]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j3VSj1j2.dat"
  histo[68]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j3VSj1j2.dat"
  histo[69]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j3VSj1j2.dat"
  histo[70]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_j3VSj1j2.dat"

  histo[71]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_j4VSj1j2.dat"
  histo[72]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j4VSj1j2.dat"
  histo[73]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j4VSj1j2.dat"
  histo[74]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j4VSj1j2.dat"
  histo[75]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_j4VSj1j2.dat"

  histo[76]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_j4VSj1j3.dat"
  histo[77]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j4VSj1j3.dat"
  histo[78]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j4VSj1j3.dat"
  histo[79]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j4VSj1j3.dat"
  histo[80]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_j4VSj1j3.dat"

  histo[81]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_j4VSj2j3.dat"
  histo[82]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_j4VSj2j3.dat"
  histo[83]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_j4VSj2j3.dat"
  histo[84]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_j4VSj2j3.dat"
  histo[85]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_j4VSj2j3.dat"

  //fillers
  histo[86]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"ystar.dat"
  histo[87]= arHistogramPtr::Create(arHistogram(-2,2,n,doSmearing));   //"ystar_normed.dat"
  histo[88]= arHistogramPtr::Create(arHistogram(-1,10,11));   //"Correlation_p1jn.dat"
  histo[89]= arHistogramPtr::Create(arHistogram(-1,10,11));   //"Correlation_p2jn.dat"
  histo[90]= arHistogramPtr::Create(arHistogram(-1,10,11));   //"Correlation_p3jn.dat"
  histo[91]= arHistogramPtr::Create(arHistogram(-1,10,11));   //"Correlation_p4jn.dat"
  histo[92]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_h.dat"
  histo[93]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_h.dat"
  histo[94]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_h.dat"
  histo[95]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_h.dat"
  histo[96]= arHistogramPtr::Create(arHistogram(0,10,10));   //"njets_inclusive.dat"
  histo[97]= arHistogramPtr::Create(arHistogram(0,10,10));   //"njets_exclusive.dat"
  histo[98]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"empty98.dat"
  histo[99]= arHistogramPtr::Create(arHistogram(-7,7,n));   //"empty99.dat"

  //corresponding partonic observables 
  histo[100]= arHistogramPtr::Create(arHistogram(0,500,n,doSmearing));   //"pt_p1.dat"
  histo[101]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p1.dat"
  histo[102]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p1.dat"
  histo[103]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p1.dat"

  histo[104]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_p2.dat"
  histo[105]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p2.dat"
  histo[106]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p2.dat"
  histo[107]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p2.dat"

  histo[108]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_p3.dat"
  histo[109]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p3.dat"
  histo[110]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p3.dat"
  histo[111]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p3.dat"

  histo[112]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_p4.dat"
  histo[113]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p4.dat"
  histo[114]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p4.dat"
  histo[115]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p4.dat"

  histo[116]= arHistogramPtr::Create(arHistogram(0,200,n,doSmearing));   //"pt_p5.dat"
  histo[117]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p5.dat"
  histo[118]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p5.dat"
  histo[119]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p5.dat"

  histo[120]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p1p2.dat"
  histo[121]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p1p2.dat"
  histo[122]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p1p2.dat"
  histo[123]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p1p2.dat"
  histo[124]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p1p2.dat"

  histo[125]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p1p3.dat"
  histo[126]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p1p3.dat"
  histo[127]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p1p3.dat"
  histo[128]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p1p3.dat"
  histo[129]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p1p3.dat"

  histo[130]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p1p4.dat"
  histo[131]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p1p4.dat"
  histo[132]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p1p4.dat"
  histo[133]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p1p4.dat"
  histo[134]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p1p4.dat"

  histo[135]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p2p3.dat"
  histo[136]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p2p3.dat"
  histo[137]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p2p3.dat"
  histo[138]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p2p3.dat"
  histo[139]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p2p3.dat"

  histo[140]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p2p4.dat"
  histo[141]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p2p4.dat"
  histo[142]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p2p4.dat"
  histo[143]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p2p4.dat"
  histo[144]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p2p3.dat"

  histo[145]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_p3p4.dat"
  histo[146]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaeta_p3p4.dat"
  histo[147]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_p3p4.dat"
  histo[148]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_p3p4.dat"
  histo[149]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_p3p4.dat"

  histo[150]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hp1.dat"
  histo[151]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hp1.dat"
  histo[152]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hp1.dat"
  histo[153]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hp1.dat"

  histo[154]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hp2.dat"
  histo[155]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hp2.dat"
  histo[156]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hp2.dat"
  histo[157]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hp2.dat"

  histo[158]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hp3.dat"
  histo[159]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hp3.dat"
  histo[160]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hp3.dat"
  histo[161]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hp3.dat"

  histo[162]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltay_hp4.dat"
  histo[163]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"deltaphi_hp4.dat"
  histo[164]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"deltaR_hp4.dat"
  histo[165]= arHistogramPtr::Create(arHistogram(0,3000,n,doSmearing));   //"m_hp4.dat"

  histo[166]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_p3VSp1p2.dat"
  histo[167]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p3VSp1p2.dat"
  histo[168]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p3VSp1p2.dat"
  histo[169]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p3VSp1p2.dat"
  histo[170]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_p3VSp1p2.dat"

  histo[171]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_p4VSp1p2.dat"
  histo[172]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p4VSp1p2.dat"
  histo[173]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p4VSp1p2.dat"
  histo[174]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p4VSp1p2.dat"
  histo[175]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_p4VSp1p2.dat"

  histo[176]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_p4VSp1p3.dat"
  histo[177]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p4VSp1p3.dat"
  histo[178]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p4VSp1p3.dat"
  histo[179]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p4VSp1p3.dat"
  histo[180]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_p4VSp1p3.dat"

  histo[181]= arHistogramPtr::Create(arHistogram(0,300,n,doSmearing));   //"pt_p4VSp2p3.dat"
  histo[182]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"y_p4VSp2p3.dat"
  histo[183]= arHistogramPtr::Create(arHistogram(-7,7,n,doSmearing));   //"eta_p4VSp2p3.dat"
  histo[184]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,doSmearing,doSmearing));   //"phi_p4VSp2p3.dat"
  histo[185]= arHistogramPtr::Create(arHistogram(0,10,n,doSmearing));   //"R_p4VSp2p3.dat"

  histo[186]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[187]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[188]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[189]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[190]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[191]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[192]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[193]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[194]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[195]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[196]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[197]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[198]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[199]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[200]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[201]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[202]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[202]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[203]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[204]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[205]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[206]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[207]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[208]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[209]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[210]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));

  histo[211]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[212]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[213]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[214]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[215]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[216]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[217]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[218]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[219]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[220]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[221]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[222]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[223]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[224]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[225]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[226]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[227]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[228]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[229]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[230]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[231]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[232]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[232]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[233]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[234]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));
  histo[235]= arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,31,doSmearing,doSmearing));

  histo2D[0]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[1]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[2]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[3]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[4]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[5]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 10, 31, 25, doSmearing,doSmearing));
  histo2D[6]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 120, 31, 24, doSmearing,doSmearing));
  histo2D[7]= ar2DHistogramPtr::Create(ar2DHistogram(-3.14159265,3.14159265, 0, 120, 31, 24, doSmearing,doSmearing));
  
  // histo[0] =arHistogramPtr::Create(arHistogram(0,10,n,true));
  // histo[1] =arHistogramPtr::Create(arHistogram(0,500,n,true));
  // histo[2] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  // histo[3] =arHistogramPtr::Create(arHistogram(0.,10.,n,true));
  // histo[4] =arHistogramPtr::Create(arHistogram(0,3000,n,true));
  // histo[5] =arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,true,true));
  // histo[6] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[7] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[8] =arHistogramPtr::Create(arHistogram(0,6,n,true));
  // histo[9] =arHistogramPtr::Create(arHistogram(0,6,n,true));
  // histo[10] =arHistogramPtr::Create(arHistogram(0,10,n,true));
  // histo[11] =arHistogramPtr::Create(arHistogram(0,500,n,true));
  // histo[12] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  // histo[13] =arHistogramPtr::Create(arHistogram(0.,10.,n,true));
  // histo[14] =arHistogramPtr::Create(arHistogram(0,3000,n,true));
  // histo[15] =arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,true,true));
  // histo[16] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[17] =arHistogramPtr::Create(arHistogram(-7,7,n,true));

  // histo[18] =arHistogramPtr::Create(arHistogram(-200,200,n,true));
  // histo[19] =arHistogramPtr::Create(arHistogram(-200,200,n,true));
  // histo[20] =arHistogramPtr::Create(arHistogram(-10,10,n,true));
  // histo[21] =arHistogramPtr::Create(arHistogram(-10,10,n,true));
  // histo[22] =arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,true,true));
  // histo[23] =arHistogramPtr::Create(arHistogram(-3.14159265,3.14159265,62,true,true));

  // histo[24] =arHistogramPtr::Create(arHistogram(-1,10,11,true));
  // histo[25] =arHistogramPtr::Create(arHistogram(-1,10,11,true));

  // histo[26] =arHistogramPtr::Create(arHistogram(1,10,n,true));
  // histo[27] =arHistogramPtr::Create(arHistogram(1,10,n,true));

  // histo[28] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[29] =arHistogramPtr::Create(arHistogram(-7,7,n,true));

  // histo[30] =arHistogramPtr::Create(arHistogram(0.,10.,n,true));

  // histo[31] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  // histo[32] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  // histo[33] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[34] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[35] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[36] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[37] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  // histo[38] =arHistogramPtr::Create(arHistogram(-2,2,n,true));
  // histo[39] =arHistogramPtr::Create(arHistogram(-100,100,n,true));



  correlarray[0][0]=0;
  correlarray[0][1]=0;
  correlarray[1][0]=0;
  correlarray[1][1]=0;
  correlarray[1][2]=0;
  Nj_delphi_25to75=0;
  Nj_delphi_0to25=0;
  Nj_delphi_75to100=0;
  Np_delphi_25to75=0;
  Np_delphi_0to25=0;
  Np_delphi_75to100=0;
}

BasicHistogramCollection::~BasicHistogramCollection(){}

void BasicHistogramCollection::writetofolder(const char * folder){
  ofstream fileoutput[236];
  ofstream fileoutput2D[8];
  ofstream gnuplotfile;
  if(strcmp(folder,"") != 0){
    mkdir(folder,0777);
    if(chdir(folder) != 0)
      cerr << "BasicHistogramCollegtion could not move into subdirectory.\n" << flush;
  }
  //Simons wishlist
  fileoutput[0].open("pt_j1.dat");
  fileoutput[1].open("eta_j1.dat");
  fileoutput[2].open("y_j1.dat");
  fileoutput[3].open("phi_j1.dat");

  fileoutput[4].open("pt_j2.dat");
  fileoutput[5].open("eta_j2.dat");
  fileoutput[6].open("y_j2.dat");
  fileoutput[7].open("phi_j2.dat");

  fileoutput[8].open("pt_j3.dat");
  fileoutput[9].open("eta_j3.dat");
  fileoutput[10].open("y_j3.dat");
  fileoutput[11].open("phi_j3.dat");

  fileoutput[12].open("pt_j4.dat");
  fileoutput[13].open("eta_j4.dat");
  fileoutput[14].open("y_j4.dat");
  fileoutput[15].open("phi_j4.dat");

  fileoutput[16].open("pt_j5.dat");
  fileoutput[17].open("eta_j5.dat");
  fileoutput[18].open("y_j5.dat");
  fileoutput[19].open("phi_j5.dat");

  fileoutput[20].open("deltay_j1j2.dat");
  fileoutput[21].open("deltaeta_j1j2.dat");
  fileoutput[22].open("deltaphi_j1j2.dat");
  fileoutput[23].open("deltaR_j1j2.dat");
  fileoutput[24].open("m_j1j2.dat");

  fileoutput[25].open("deltay_j1j3.dat");
  fileoutput[26].open("deltaeta_j1j3.dat");
  fileoutput[27].open("deltaphi_j1j3.dat");
  fileoutput[28].open("deltaR_j1j3.dat");
  fileoutput[29].open("m_j1j3.dat");

  fileoutput[30].open("deltay_j1j4.dat");
  fileoutput[31].open("deltaeta_j1j4.dat");
  fileoutput[32].open("deltaphi_j1j4.dat");
  fileoutput[33].open("deltaR_j1j4.dat");
  fileoutput[34].open("m_j1j4.dat");

  fileoutput[35].open("deltay_j2j3.dat");
  fileoutput[36].open("deltaeta_j2j3.dat");
  fileoutput[37].open("deltaphi_j2j3.dat");
  fileoutput[38].open("deltaR_j2j3.dat");
  fileoutput[39].open("m_j2j3.dat");

  fileoutput[40].open("deltay_j2j4.dat");
  fileoutput[41].open("deltaeta_j2j4.dat");
  fileoutput[42].open("deltaphi_j2j4.dat");
  fileoutput[43].open("deltaR_j2j4.dat");
  fileoutput[44].open("m_j2j4.dat");

  fileoutput[45].open("deltay_j3j4.dat");
  fileoutput[46].open("deltaeta_j3j4.dat");
  fileoutput[47].open("deltaphi_j3j4.dat");
  fileoutput[48].open("deltaR_j3j4.dat");
  fileoutput[49].open("m_j3j4.dat");

  fileoutput[50].open("deltay_hj1.dat");
  fileoutput[51].open("deltaphi_hj1.dat");
  fileoutput[52].open("deltaR_hj1.dat");
  fileoutput[53].open("m_hj1.dat");

  fileoutput[54].open("deltay_hj2.dat");
  fileoutput[55].open("deltaphi_hj2.dat");
  fileoutput[56].open("deltaR_hj2.dat");
  fileoutput[57].open("m_hj2.dat");

  fileoutput[58].open("deltay_hj3.dat");
  fileoutput[59].open("deltaphi_hj3.dat");
  fileoutput[60].open("deltaR_hj3.dat");
  fileoutput[61].open("m_hj3.dat");

  fileoutput[62].open("deltay_hj4.dat");
  fileoutput[63].open("deltaphi_hj4.dat");
  fileoutput[64].open("deltaR_hj4.dat");
  fileoutput[65].open("m_hj4.dat");

  fileoutput[66].open("pt_j3VSj1j2.dat");
  fileoutput[67].open("y_j3VSj1j2.dat");
  fileoutput[68].open("eta_j3VSj1j2.dat");
  fileoutput[69].open("phi_j3VSj1j2.dat");
  fileoutput[70].open("R_j3VSj1j2.dat");

  fileoutput[71].open("pt_j4VSj1j2.dat");
  fileoutput[72].open("y_j4VSj1j2.dat");
  fileoutput[73].open("eta_j4VSj1j2.dat");
  fileoutput[74].open("phi_j4VSj1j2.dat");
  fileoutput[75].open("R_j4VSj1j2.dat");

  fileoutput[76].open("pt_j4VSj1j3.dat");
  fileoutput[77].open("y_j4VSj1j3.dat");
  fileoutput[78].open("eta_j4VSj1j3.dat");
  fileoutput[79].open("phi_j4VSj1j3.dat");
  fileoutput[80].open("R_j4VSj1j3.dat");

  fileoutput[81].open("pt_j4VSj2j3.dat");
  fileoutput[82].open("y_j4VSj2j3.dat");
  fileoutput[83].open("eta_j4VSj2j3.dat");
  fileoutput[84].open("phi_j4VSj2j3.dat");
  fileoutput[85].open("R_j4VSj2j3.dat");

  //fillers
  fileoutput[86].open("ystar.dat");
  fileoutput[87].open("ystar_normed.dat");
  fileoutput[88].open("Correlation_p1jn.dat");
  fileoutput[89].open("Correlation_p2jn.dat");
  fileoutput[90].open("Correlation_p3jn.dat");
  fileoutput[91].open("Correlation_p4jn.dat");
  fileoutput[92].open("pt_h.dat");
  fileoutput[93].open("eta_h.dat");
  fileoutput[94].open("y_h.dat");
  fileoutput[95].open("phi_h.dat");
  fileoutput[96].open("njets_inclusive.dat");
  fileoutput[97].open("njets_exclusive.dat");
  fileoutput[98].open("nevents_eta3_smeared.dat");
  fileoutput[99].open("nevents_eta3.dat");

  //corresponding partonic observables 
  fileoutput[100].open("pt_p1.dat");
  fileoutput[101].open("eta_p1.dat");
  fileoutput[102].open("y_p1.dat");
  fileoutput[103].open("phi_p1.dat");

  fileoutput[104].open("pt_p2.dat");
  fileoutput[105].open("eta_p2.dat");
  fileoutput[106].open("y_p2.dat");
  fileoutput[107].open("phi_p2.dat");

  fileoutput[108].open("pt_p3.dat");
  fileoutput[109].open("eta_p3.dat");
  fileoutput[110].open("y_p3.dat");
  fileoutput[111].open("phi_p3.dat");

  fileoutput[112].open("pt_p4.dat");
  fileoutput[113].open("eta_p4.dat");
  fileoutput[114].open("y_p4.dat");
  fileoutput[115].open("phi_p4.dat");

  fileoutput[116].open("pt_p5.dat");
  fileoutput[117].open("eta_p5.dat");
  fileoutput[118].open("y_p5.dat");
  fileoutput[119].open("phi_p5.dat");

  fileoutput[120].open("deltay_p1p2.dat");
  fileoutput[121].open("deltaeta_p1p2.dat");
  fileoutput[122].open("deltaphi_p1p2.dat");
  fileoutput[123].open("deltaR_p1p2.dat");
  fileoutput[124].open("m_p1p2.dat");

  fileoutput[125].open("deltay_p1p3.dat");
  fileoutput[126].open("deltaeta_p1p3.dat");
  fileoutput[127].open("deltaphi_p1p3.dat");
  fileoutput[128].open("deltaR_p1p3.dat");
  fileoutput[129].open("m_p1p3.dat");

  fileoutput[130].open("deltay_p1p4.dat");
  fileoutput[131].open("deltaeta_p1p4.dat");
  fileoutput[132].open("deltaphi_p1p4.dat");
  fileoutput[133].open("deltaR_p1p4.dat");
  fileoutput[134].open("m_p1p4.dat");

  fileoutput[135].open("deltay_p2p3.dat");
  fileoutput[136].open("deltaeta_p2p3.dat");
  fileoutput[137].open("deltaphi_p2p3.dat");
  fileoutput[138].open("deltaR_p2p3.dat");
  fileoutput[139].open("m_p2p3.dat");

  fileoutput[140].open("deltay_p2p4.dat");
  fileoutput[141].open("deltaeta_p2p4.dat");
  fileoutput[142].open("deltaphi_p2p4.dat");
  fileoutput[143].open("deltaR_p2p4.dat");
  fileoutput[144].open("m_p2p4.dat");

  fileoutput[145].open("deltay_p3p4.dat");
  fileoutput[146].open("deltaeta_p3p4.dat");
  fileoutput[147].open("deltaphi_p3p4.dat");
  fileoutput[148].open("deltaR_p3p4.dat");
  fileoutput[149].open("m_p3p4.dat");

  fileoutput[150].open("deltay_hp1.dat");
  fileoutput[151].open("deltaphi_hp1.dat");
  fileoutput[152].open("deltaR_hp1.dat");
  fileoutput[153].open("m_hp1.dat");

  fileoutput[154].open("deltay_hp2.dat");
  fileoutput[155].open("deltaphi_hp2.dat");
  fileoutput[156].open("deltaR_hp2.dat");
  fileoutput[157].open("m_hp2.dat");

  fileoutput[158].open("deltay_hp3.dat");
  fileoutput[159].open("deltaphi_hp3.dat");
  fileoutput[160].open("deltaR_hp3.dat");
  fileoutput[161].open("m_hp3.dat");

  fileoutput[162].open("deltay_hp4.dat");
  fileoutput[163].open("deltaphi_hp4.dat");
  fileoutput[164].open("deltaR_hp4.dat");
  fileoutput[165].open("m_hp4.dat");

  fileoutput[166].open("pt_p3VSp1p2.dat");
  fileoutput[167].open("y_p3VSp1p2.dat");
  fileoutput[168].open("eta_p3VSp1p2.dat");
  fileoutput[169].open("phi_p3VSp1p2.dat");
  fileoutput[170].open("R_p3VSp1p2.dat");

  fileoutput[171].open("pt_p4VSp1p2.dat");
  fileoutput[172].open("y_p4VSp1p2.dat");
  fileoutput[173].open("eta_p4VSp1p2.dat");
  fileoutput[174].open("phi_p4VSp1p2.dat");
  fileoutput[175].open("R_p4VSp1p2.dat");

  fileoutput[176].open("pt_p4VSp1p3.dat");
  fileoutput[177].open("y_p4VSp1p3.dat");
  fileoutput[178].open("eta_p4VSp1p3.dat");
  fileoutput[179].open("phi_p4VSp1p3.dat");
  fileoutput[180].open("R_p4VSp1p3.dat");

  fileoutput[181].open("pt_p4VSp2p3.dat");
  fileoutput[182].open("y_p4VSp2p3.dat");
  fileoutput[183].open("eta_p4VSp2p3.dat");
  fileoutput[184].open("phi_p4VSp2p3.dat");
  fileoutput[185].open("R_p4VSp2p3.dat");

  fileoutput[186].open("phi_j1j2eta_j1j2.004.2D.dat");
  fileoutput[187].open("phi_j1j2eta_j1j2.008.2D.dat");
  fileoutput[188].open("phi_j1j2eta_j1j2.012.2D.dat");
  fileoutput[189].open("phi_j1j2eta_j1j2.016.2D.dat");
  fileoutput[190].open("phi_j1j2eta_j1j2.020.2D.dat");
  fileoutput[191].open("phi_j1j2eta_j1j2.024.2D.dat");
  fileoutput[192].open("phi_j1j2eta_j1j2.028.2D.dat");
  fileoutput[193].open("phi_j1j2eta_j1j2.032.2D.dat");
  fileoutput[194].open("phi_j1j2eta_j1j2.036.2D.dat");
  fileoutput[195].open("phi_j1j2eta_j1j2.040.2D.dat");
  fileoutput[196].open("phi_j1j2eta_j1j2.044.2D.dat");
  fileoutput[197].open("phi_j1j2eta_j1j2.048.2D.dat");
  fileoutput[198].open("phi_j1j2eta_j1j2.052.2D.dat");
  fileoutput[199].open("phi_j1j2eta_j1j2.056.2D.dat");
  fileoutput[200].open("phi_j1j2eta_j1j2.060.2D.dat");
  fileoutput[201].open("phi_j1j2eta_j1j2.064.2D.dat");
  fileoutput[202].open("phi_j1j2eta_j1j2.068.2D.dat");
  fileoutput[203].open("phi_j1j2eta_j1j2.072.2D.dat");
  fileoutput[204].open("phi_j1j2eta_j1j2.076.2D.dat");
  fileoutput[205].open("phi_j1j2eta_j1j2.080.2D.dat");
  fileoutput[206].open("phi_j1j2eta_j1j2.084.2D.dat");
  fileoutput[207].open("phi_j1j2eta_j1j2.088.2D.dat");
  fileoutput[208].open("phi_j1j2eta_j1j2.092.2D.dat");
  fileoutput[209].open("phi_j1j2eta_j1j2.096.2D.dat");
  fileoutput[210].open("phi_j1j2eta_j1j2.100.2D.dat");

  fileoutput[211].open("phi_j2j3eta_j2j3.004.2D.dat");
  fileoutput[212].open("phi_j2j3eta_j2j3.008.2D.dat");
  fileoutput[213].open("phi_j2j3eta_j2j3.012.2D.dat");
  fileoutput[214].open("phi_j2j3eta_j2j3.016.2D.dat");
  fileoutput[215].open("phi_j2j3eta_j2j3.020.2D.dat");
  fileoutput[216].open("phi_j2j3eta_j2j3.024.2D.dat");
  fileoutput[217].open("phi_j2j3eta_j2j3.028.2D.dat");
  fileoutput[218].open("phi_j2j3eta_j2j3.032.2D.dat");
  fileoutput[219].open("phi_j2j3eta_j2j3.036.2D.dat");
  fileoutput[220].open("phi_j2j3eta_j2j3.040.2D.dat");
  fileoutput[221].open("phi_j2j3eta_j2j3.044.2D.dat");
  fileoutput[222].open("phi_j2j3eta_j2j3.048.2D.dat");
  fileoutput[223].open("phi_j2j3eta_j2j3.052.2D.dat");
  fileoutput[224].open("phi_j2j3eta_j2j3.056.2D.dat");
  fileoutput[225].open("phi_j2j3eta_j2j3.060.2D.dat");
  fileoutput[226].open("phi_j2j3eta_j2j3.064.2D.dat");
  fileoutput[227].open("phi_j2j3eta_j2j3.068.2D.dat");
  fileoutput[228].open("phi_j2j3eta_j2j3.072.2D.dat");
  fileoutput[229].open("phi_j2j3eta_j2j3.076.2D.dat");
  fileoutput[230].open("phi_j2j3eta_j2j3.080.2D.dat");
  fileoutput[231].open("phi_j2j3eta_j2j3.084.2D.dat");
  fileoutput[232].open("phi_j2j3eta_j2j3.088.2D.dat");
  fileoutput[233].open("phi_j2j3eta_j2j3.092.2D.dat");
  fileoutput[234].open("phi_j2j3eta_j2j3.096.2D.dat");
  fileoutput[235].open("phi_j2j3eta_j2j3.100.2D.dat");

  fileoutput2D[0].open("phi_j1j2eta_j1j2.dat.2D");
  fileoutput2D[1].open("phi_j2j3eta_j2j3.dat.2D");
  fileoutput2D[2].open("phi_j1jSameeta_j1jSame.dat.2D");
  fileoutput2D[3].open("phi_j1jOppositeeta_j1jOpposite.dat.2D");
  fileoutput2D[4].open("phi_j2jSameeta_j2jSame.dat.2D");
  fileoutput2D[5].open("phi_j2jOppositeeta_j2jOpposite.dat.2D");
  fileoutput2D[6].open("phi_j1j2pt_j3.dat.2D");
  fileoutput2D[7].open("phi_j2j3pt_j3.dat.2D");

 // fileoutput[0].open("R_jj");
  // fileoutput[1].open("pt_j1");
  // fileoutput[2].open("pt_j2");
  // fileoutput[3].open("deltaeta_jj");
  // fileoutput[4].open("m_jj");
  // fileoutput[5].open("phi_jj");
  // fileoutput[6].open("eta_j1");
  // fileoutput[7].open("eta_j2");

  // fileoutput[8].open("R_pj1");
  // fileoutput[9].open("R_pj2");

  // fileoutput[10].open("R_pp");
  // fileoutput[11].open("pt_p1");
  // fileoutput[12].open("pt_p2");
  // fileoutput[13].open("deltay_pp");
  // fileoutput[14].open("m_pp");
  // fileoutput[15].open("phi_pp");
  // fileoutput[16].open("y_p1");
  // fileoutput[17].open("y_p2");

  // fileoutput[18].open("deltapt_pj1");
  // fileoutput[19].open("deltapt_pj2");
  // fileoutput[20].open("deltay_pj1");
  // fileoutput[21].open("deltay_pj2");
  // fileoutput[22].open("deltaphi_pj1");
  // fileoutput[23].open("deltaphi_pj2");

  // fileoutput[24].open("n_p1jnCorrelation");
  // fileoutput[25].open("n_p2jnCorrelation");

  // fileoutput[26].open("pt_j1:pt_j2");
  // fileoutput[27].open("pt_p1:pt_p2");

  // fileoutput[28].open("y_j1");
  // fileoutput[29].open("y_j2");

  // fileoutput[30].open("deltay_jj");

  
  // fileoutput[31].open("pt_j3");
  // fileoutput[32].open("pt_j4");
  // fileoutput[33].open("eta_j3");
  // fileoutput[34].open("eta_j4");
  // fileoutput[35].open("y_j3");
  // fileoutput[36].open("y_j4");
  // fileoutput[37].open("y_star");
  // fileoutput[38].open("y_star_normed");
  // fileoutput[39].open("weights");

  for (int i=0; i<=235; i++){
    (*histo[i]).simpleOutput(fileoutput[i],true,true,true);
    fileoutput[i].close();
  }
  for (int i=0; i<=7; i++){
    (*histo2D[i]).simpleOutput(fileoutput2D[i],true,true,true);
    fileoutput2D[i].close();
  }

  ofstream Aphiout;
  Aphiout.open("A_phi");
  Aphiout<<"A_phi_jets   =   "<<A_phi_jets<<endl
  	 <<"A_phi_partons=   "<<A_phi_partons; 
  
  ofstream correlationout;
  correlationout.open("correlsample.gp");
  
  correlationout<<"set output \""<<folder<<".ps\""<<endl
  		<<"set terminal postscript enhanced color"<<endl
  		<<"#set style data boxes"<<endl
  		<<"#set style data histeps"<<endl
  		<<"set style data fsteps"<<endl
  		<<"set key right top"<<endl
  		<<endl
  		<<"set title \""<<folder<<" Correlation samples\""<<endl
  		<<"set xrange [0.0:5.0]"<<endl
  		<<"set ylabel \"# events\""<<endl
  		<<"set xlabel \" sample\""<<endl
  		<<"set xtics (\"j1=p1 & j2=p2\" 0.5, \"j1=p2 & j2=p1\" 1.5, \"j1=p1 & j2!=p2\" 2.5, \"j1!=p1 & j2=p2\" 3.5, \"j1!=p1/2 & j2=p1/2\" 4.5)    "<<endl
  		<<"plot \"-\""<<endl
  		<<"0   0"<<endl
  		<<"1   "<<correlarray[0][0]<<endl
  		<<"2   "<<correlarray[0][2]<<endl
  		<<"3   "<<correlarray[0][1]<<endl
  		<<"4   "<<correlarray[1][0]<<endl
  		<<"5   "<<correlarray[1][1]<<endl
  		<<endl;


  gnuplotfile.open("BasicHistograms.gp");
  gnuplotfile<<"set output \""<<folder<<".ps\""<<endl
  	     <<"set terminal postscript enhanced color"<<endl
  	     <<"#set style data boxes"<<endl
  	     <<"#set style data histeps"<<endl
  	     <<"set style data steps"<<endl
  	     <<"set key right top"<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of hardest tagging jet\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,j_1} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,j_1} [GeV]\""<<endl
  	     <<"plot 	\"pt_j1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of hardest tagging jet {/Symbol h}_{j_1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_1} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{j_1}\""<<endl
  	     <<"plot 	\"eta_j1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of hardest tagging jet {y}_{j_1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{j_1} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{j_1}\""<<endl
  	     <<"plot 	\"y_j1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of hardest tagging jet {y}_{j_1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j_1} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j_1}\""<<endl
  	     <<"plot 	\"phi_j1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of 2nd hardest tagging jet\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,j_2} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,j_2} [GeV]\""<<endl
  	     <<"plot 	\"pt_j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of 2nd hardest tagging jet {/Symbol h}_{j_2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_2} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{j_2}\""<<endl
  	     <<"plot 	\"eta_j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of 2nd hardest tagging jet {y}_{j_2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{j_2} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{j_2}\""<<endl
  	     <<"plot 	\"y_j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of 2nd hardest tagging jet {y}_{j_2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j_2} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j_2}\""<<endl
  	     <<"plot 	\"phi_j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of 3nd hardest tagging jet\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,j_3} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,j_3} [GeV]\""<<endl
  	     <<"plot 	\"pt_j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of 3nd hardest tagging jet {/Symbol h}_{j_3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{j_3}\""<<endl
  	     <<"plot 	\"eta_j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of 3nd hardest tagging jet {y}_{j_3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{j_3} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{j_3}\""<<endl
  	     <<"plot 	\"y_j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of 3nd hardest tagging jet {y}_{j_3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j_3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j_3}\""<<endl
  	     <<"plot 	\"phi_j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of 4nd hardest tagging jet\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,j_4} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,j_4} [GeV]\""<<endl
  	     <<"plot 	\"pt_j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of 4nd hardest tagging jet {/Symbol h}_{j_4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{j_4}\""<<endl
  	     <<"plot 	\"eta_j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of 4nd hardest tagging jet {y}_{j_4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{j_4} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{j_4}\""<<endl
  	     <<"plot 	\"y_j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of 4nd hardest tagging jet {y}_{j_4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j_4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j_4}\""<<endl
  	     <<"plot 	\"phi_j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of 5nd hardest tagging jet\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,j_5} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,j_5} [GeV]\""<<endl
  	     <<"plot 	\"pt_j5.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p5.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of 5nd hardest tagging jet {/Symbol h}_{j_5}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_5} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{j_5}\""<<endl
  	     <<"plot 	\"eta_j5.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p5.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of 5nd hardest tagging jet {y}_{j_5}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{j_5} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{j_5}\""<<endl
  	     <<"plot 	\"y_j5.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p5.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of 5nd hardest tagging jet {y}_{j_5}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j_5} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j_5}\""<<endl
  	     <<"plot 	\"phi_j5.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p5.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j1j2}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j1j2} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j1j2}\""<<endl
             <<"plot    \"deltay_j1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p1p2.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j1j2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j1j2} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j1j2}\""<<endl
  	     <<"plot 	\"deltaeta_j1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j1j2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j1j2} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j1j2}\""<<endl
  	     <<"plot 	\"deltaphi_j1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j1j2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j1j2} [fb]\""<<endl
  	     <<"set xlabel \"R_{j1j2}\""<<endl
  	     <<"plot 	\"deltaR_j1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j1j2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j1j2} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j1j2} [GeV]\""<<endl
  	     <<"plot 	\"m_j1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j1j3}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j1j3} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j1j3}\""<<endl
             <<"plot    \"deltay_j1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p1p3.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j1j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j1j3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j1j3}\""<<endl
  	     <<"plot 	\"deltaeta_j1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j1j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j1j3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j1j3}\""<<endl
  	     <<"plot 	\"deltaphi_j1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j1j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j1j3} [fb]\""<<endl
  	     <<"set xlabel \"R_{j1j3}\""<<endl
  	     <<"plot 	\"deltaR_j1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j1j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j1j3} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j1j3} [GeV]\""<<endl
  	     <<"plot 	\"m_j1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j1j4}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j1j4} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j1j4}\""<<endl
             <<"plot    \"deltay_j1j4.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p1p4.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j1j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j1j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j1j4}\""<<endl
  	     <<"plot 	\"deltaeta_j1j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p1p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j1j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j1j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j1j4}\""<<endl
  	     <<"plot 	\"deltaphi_j1j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p1p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j1j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j1j4} [fb]\""<<endl
  	     <<"set xlabel \"R_{j1j4}\""<<endl
  	     <<"plot 	\"deltaR_j1j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p1p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j1j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j1j4} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j1j4} [GeV]\""<<endl
  	     <<"plot 	\"m_j1j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p1p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j2j3}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j2j3} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j2j3}\""<<endl
             <<"plot    \"deltay_j2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p2p3.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j2j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j2j3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j2j3}\""<<endl
  	     <<"plot 	\"deltaeta_j2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j2j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j2j3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j2j3}\""<<endl
  	     <<"plot 	\"deltaphi_j2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j2j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j2j3} [fb]\""<<endl
  	     <<"set xlabel \"R_{j2j3}\""<<endl
  	     <<"plot 	\"deltaR_j2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j2j3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j2j3} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j2j3} [GeV]\""<<endl
  	     <<"plot 	\"m_j2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j2j4}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j2j4} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j2j4}\""<<endl
             <<"plot    \"deltay_j2j4.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p2p4.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j2j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j2j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j2j4}\""<<endl
  	     <<"plot 	\"deltaeta_j2j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p2p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j2j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j2j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j2j4}\""<<endl
  	     <<"plot 	\"deltaphi_j2j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p2p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j2j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j2j4} [fb]\""<<endl
  	     <<"set xlabel \"R_{j2j4}\""<<endl
  	     <<"plot 	\"deltaR_j2j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p2p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j2j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j2j4} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j2j4} [GeV]\""<<endl
  	     <<"plot 	\"m_j2j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p2p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{j3j4}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{j3j4} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{j3j4}\""<<endl
             <<"plot    \"deltay_j3j4.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_p3p4.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{j3j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{j3j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol D}{/Symbol h}_{j3j4}\""<<endl
  	     <<"plot 	\"deltaeta_j3j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaeta_p3p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{j3j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{j3j4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{j3j4}\""<<endl
  	     <<"plot 	\"deltaphi_j3j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_p3p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{j3j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{j3j4} [fb]\""<<endl
  	     <<"set xlabel \"R_{j3j4}\""<<endl
  	     <<"plot 	\"deltaR_j3j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_p3p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{j3j4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{j3j4} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{j3j4} [GeV]\""<<endl
  	     <<"plot 	\"m_j3j4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_p3p4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{hj1}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{hj1} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{hj1}\""<<endl
             <<"plot    \"deltay_hj1.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_hp1.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{hj1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{hj1} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{hj1}\""<<endl
  	     <<"plot 	\"deltaphi_hj1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_hp1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{hj1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{hj1} [fb]\""<<endl
  	     <<"set xlabel \"R_{hj1}\""<<endl
  	     <<"plot 	\"deltaR_hj1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_hp1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{hj1}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{hj1} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{hj1} [GeV]\""<<endl
  	     <<"plot 	\"m_hj1.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_hp1.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{hj2}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{hj2} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{hj2}\""<<endl
             <<"plot    \"deltay_hj2.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_hp2.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{hj2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{hj2} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{hj2}\""<<endl
  	     <<"plot 	\"deltaphi_hj2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_hp2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{hj2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{hj2} [fb]\""<<endl
  	     <<"set xlabel \"R_{hj2}\""<<endl
  	     <<"plot 	\"deltaR_hj2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_hp2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{hj2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{hj2} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{hj2} [GeV]\""<<endl
  	     <<"plot 	\"m_hj2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_hp2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{hj3}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{hj3} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{hj3}\""<<endl
             <<"plot    \"deltay_hj3.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_hp3.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{hj3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{hj3} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{hj3}\""<<endl
  	     <<"plot 	\"deltaphi_hj3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_hp3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{hj3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{hj3} [fb]\""<<endl
  	     <<"set xlabel \"R_{hj3}\""<<endl
  	     <<"plot 	\"deltaR_hj3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_hp3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{hj3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{hj3} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{hj3} [GeV]\""<<endl
  	     <<"plot 	\"m_hj3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_hp3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
             <<"set title \""<<folder<<": {/Symbol D}y_{hj4}\""<<endl
             <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{hj4} [fb]\""<<endl
             <<"set xlabel \"{/Symbol D}y_{hj4}\""<<endl
             <<"plot    \"deltay_hj4.dat\" using 1:7 title \"jet level\",\\"<<endl
             <<"        \"deltay_hp4.dat\" using 1:7 title \"associated partons\""<<endl
             <<endl
  	     <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{hj4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{hj4} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{hj4}\""<<endl
  	     <<"plot 	\"deltaphi_hj4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaphi_hp4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": legoplot separation R_{hj4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{hj4} [fb]\""<<endl
  	     <<"set xlabel \"R_{hj4}\""<<endl
  	     <<"plot 	\"deltaR_hj4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"deltaR_hp4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": dijet invariant mass m_{hj4}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dm_{hj4} [fb/GeV]\""<<endl
  	     <<"set xlabel \"m_{hj4} [GeV]\""<<endl
  	     <<"plot 	\"m_hj4.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"m_hp4.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_{T,j_3} - (p_{T,j_1} + p_{T,j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,3VS12} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,3VS12} [GeV]\""<<endl
  	     <<"plot 	\"pt_j3VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p3VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol h}_{j_3} - ({/Symbol h}_{j_1} + {/Symbol h}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{3VS12} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{3VS12}\""<<endl
  	     <<"plot 	\"eta_j3VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p3VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {y}_{j_3} - ({y}_{j_1} + {y}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{3VS12} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{3VS12}\""<<endl
  	     <<"plot 	\"y_j3VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p3VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol F}_{j_3} - ({/Symbol F}_{j_1} + {/Symbol F}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{3VS12} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{3VS12}\""<<endl
  	     <<"plot 	\"phi_j3VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p3VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": R_{3,1+2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{3VS12} [fb]\""<<endl
  	     <<"set xlabel \"R_{3VS12}\""<<endl
  	     <<"plot 	\"R_j3VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"R_p3VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_{T,j_4} - (p_{T,j_1} + p_{T,j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,4VS12} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,4VS12} [GeV]\""<<endl
  	     <<"plot 	\"pt_j4VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p4VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol h}_{j_4} - ({/Symbol h}_{j_1} + {/Symbol h}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{4VS12} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{4VS12}\""<<endl
  	     <<"plot 	\"eta_j4VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p4VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {y}_{j_4} - ({y}_{j_1} + {y}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{4VS12} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{4VS12}\""<<endl
  	     <<"plot 	\"y_j4VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p4VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol F}_{j_4} - ({/Symbol F}_{j_1} + {/Symbol F}_{j_2})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{4VS12} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{4VS12}\""<<endl
  	     <<"plot 	\"phi_j4VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p4VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": R_{4,1+2}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{4VS12} [fb]\""<<endl
  	     <<"set xlabel \"R_{4VS12}\""<<endl
  	     <<"plot 	\"R_j4VSj1j2.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"R_p4VSp1p2.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_{T,j_4} - (p_{T,j_1} + p_{T,j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,4VS13} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,4VS13} [GeV]\""<<endl
  	     <<"plot 	\"pt_j4VSj1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p4VSp1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol h}_{j_4} - ({/Symbol h}_{j_1} + {/Symbol h}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{4VS13} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{4VS13}\""<<endl
  	     <<"plot 	\"eta_j4VSj1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p4VSp1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {y}_{j_4} - ({y}_{j_1} + {y}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{4VS13} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{4VS13}\""<<endl
  	     <<"plot 	\"y_j4VSj1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p4VSp1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol F}_{j_4} - ({/Symbol F}_{j_1} + {/Symbol F}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{4VS13} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{4VS13}\""<<endl
  	     <<"plot 	\"phi_j4VSj1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p4VSp1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": R_{4,1+3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{4VS13} [fb]\""<<endl
  	     <<"set xlabel \"R_{4VS13}\""<<endl
  	     <<"plot 	\"R_j4VSj1j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"R_p4VSp1p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_{T,j_4} - (p_{T,j_2} + p_{T,j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,4VS23} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,4VS23} [GeV]\""<<endl
  	     <<"plot 	\"pt_j4VSj2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"pt_p4VSp2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol h}_{j_4} - ({/Symbol h}_{j_2} + {/Symbol h}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{4VS23} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{4VS23}\""<<endl
  	     <<"plot 	\"eta_j4VSj2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"eta_p4VSp2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {y}_{j_4} - ({y}_{j_2} + {y}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{4VS23} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{4VS23}\""<<endl
  	     <<"plot 	\"y_j4VSj2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"y_p4VSp2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": {/Symbol F}_{j_4} - ({/Symbol F}_{j_2} + {/Symbol F}_{j_3})\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{4VS23} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{4VS23}\""<<endl
  	     <<"plot 	\"phi_j4VSj2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"phi_p4VSp2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": R_{4,2+3}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dR_{4VS23} [fb]\""<<endl
  	     <<"set xlabel \"R_{4VS23}\""<<endl
  	     <<"plot 	\"R_j4VSj2j3.dat\" using 1:7 title \"jet level\",\\"<<endl
  	     <<"	\"R_p4VSp2p3.dat\" using 1:7 title \"associated partons\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": y^*\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y^*} [fb]\""<<endl
  	     <<"set xlabel \"y^*\""<<endl
  	     <<"plot 	\"ystar.dat\" using 1:7 title \"jet level\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": ||y^*||\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{||y^*||} [fb]\""<<endl
  	     <<"set xlabel \"||y^*||\""<<endl
  	     <<"plot 	\"ystar_normed.dat\" using 1:7 title \"jet level\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": p_T of  higgs\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dp_{T,h} [fb/GeV]\""<<endl
  	     <<"set xlabel \"p_{T,h} [GeV]\""<<endl
  	     <<"plot 	\"pt_h.dat\" using 1:7 title \"higgs\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": pseudorapidity of higgs {/Symbol h}_{h}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{h} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol h}_{h}\""<<endl
  	     <<"plot 	\"eta_h.dat\" using 1:7 title \"higgs\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": rapidity of higgs {y}_{h}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{y}_{h} [fb]\""<<endl
  	     <<"set xlabel \"{y}_{h}\""<<endl
  	     <<"plot 	\"y_h.dat\" using 1:7 title \"higgs\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": azimuthal angle of higgs {y}_{h}\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{h} [fb]\""<<endl
  	     <<"set xlabel \"{/Symbol F}_{h}\""<<endl
  	     <<"plot 	\"phi_h.dat\" using 1:7 title \"higgs\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": inclusive number of jets\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dn_{inclusive} [fb]\""<<endl
  	     <<"set xlabel \"n_{inclusive}\""<<endl
  	     <<"plot 	\"njets_inclusive.dat\" using 2:7 title \"inclusive\""<<endl
  	     <<endl
  	     <<"set title \""<<folder<<": exclusive number of jets\""<<endl
  	     <<"set ylabel \"d{/Symbol s}/dn_{exclusive} [fb]\""<<endl
  	     <<"set xlabel \"n_{exclusive}\""<<endl
  	     <<"plot 	\"njets_exclusive.dat\" using 2:7 title \"exclusive\""<<endl
  	     <<endl
  	     <<"set style data fsteps"<<endl
  	     <<"set title \""<<folder<<": Correlation plot\""<<endl
  	     <<"set ylabel \"relative frequency\""<<endl
  	     <<"set xlabel \"jet #\""<<endl
  	     <<"set xrange [-1.5:8.5]"<<endl
  	     <<"plot 	\"Correlation_p1jn.dat\" using 2:4 title \"correlated jet of hardest parton\",\\"<<endl
  	     <<"	\"Correlation_p2jn.dat\" using 2:4 title \"correlated jet of 2nd hardest parton\",\\"<<endl
	     <<"	\"Correlation_p3jn.dat\" using 2:4 title \"correlated jet of 2nd hardest parton\",\\"<<endl
	     <<"	\"Correlation_p4jn.dat\" using 2:4 title \"correlated jet of 2nd hardest parton\""<<endl;






  	     // <<"set title \""<<folder<<": azimuthal angle distribution {/Symbol F}_{jj}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol F}_{jj} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol F}_{jj}\""<<endl
  	     // <<"plot 	\"phi_jj\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"phi_pp\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{jj}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{jj} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}{/Symbol h}_{jj}\""<<endl
  	     // <<"plot 	\"deltaeta_jj\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"deltay_pp\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
             // <<"set title \""<<folder<<": {/Symbol D}y_{jj}\""<<endl
             // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}y_{jj} [fb]\""<<endl
             // <<"set xlabel \"{/Symbol D}y_{jj}\""<<endl
             // <<"plot    \"deltay_jj\" using 1:7 title \"jet level\",\\"<<endl
             // <<"        \"deltay_pp\" using 1:7 title \"associated partons\""<<endl
             // <<endl
  	     // <<"set title \""<<folder<<": y^*\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{y^*} [fb]\""<<endl
  	     // <<"set xlabel \"y^*\""<<endl
  	     // <<"plot 	\"y_star\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": ||y^*||\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{||y^*||} [fb]\""<<endl
  	     // <<"set xlabel \"||y^*||\""<<endl
  	     // <<"plot 	\"y_star_normed\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": pseudorapidity of 2nd hardest tagging jet {/Symbol h}_{j_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_2} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol h}_{j_2}\""<<endl
  	     // <<"plot 	\"eta_j2\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"y_p2\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": pseudorapidity of 3rd hardest tagging jet {/Symbol h}_{j_3}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_3} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol h}_{j_3}\""<<endl
  	     // <<"plot 	\"eta_j3\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": pseudorapidity of 4th hardest tagging jet {/Symbol h}_{j_4}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j_4} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol h}_{j_4}\""<<endl
  	     // <<"plot 	\"eta_j4\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": rapidity of 2nd hardest tagging jet {y}_{j_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{y}_{j_2} [fb]\""<<endl
  	     // <<"set xlabel \"{y}_{j_2}\""<<endl
  	     // <<"plot 	\"y_j2\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"y_p2\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": rapidity of 3rd hardest tagging jet {y}_{j_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{y}_{j_3} [fb]\""<<endl
  	     // <<"set xlabel \"{y}_{j_3}\""<<endl
  	     // <<"plot 	\"y_j3\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": rapidity of 4th hardest tagging jet {y}_{j_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{y}_{j_4} [fb]\""<<endl
  	     // <<"set xlabel \"{y}_{j_4}\""<<endl
  	     // <<"plot 	\"y_j4\" using 1:7 title \"jet level\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": dijet invariant mass m_{jj}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/dm_{jj} [fb/MeV]\""<<endl
  	     // <<"set xlabel \"m_{jj} [MeV]\""<<endl
  	     // <<"plot 	\"m_jj\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"m_pp\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}p_{T,pj_1}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}p_{T,pj_1} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}p_{T,pj_1}\""<<endl
  	     // <<"plot 	\"deltapt_pj1\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}p_{T,pj_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}p_{T,pj_2} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}p_{T,pj_2}\""<<endl
  	     // <<"plot 	\"deltapt_pj2\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{pj_1}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{pj_1} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}{/Symbol h}_{pj_1}\""<<endl
  	     // <<"plot 	\"deltay_pj1\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}{/Symbol h}_{pj_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol h}_{pj_2} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}{/Symbol h}_{pj_2}\""<<endl
  	     // <<"plot 	\"deltay_pj2\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{pj_1}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol F}_{pj_1} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}{/Symbol F}_{pj_1}\""<<endl
  	     // <<"plot 	\"deltaphi_pj1\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": {/Symbol D}{/Symbol F}_{pj_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d{/Symbol D}{/Symbol F}_{pj_2} [fb]\""<<endl
  	     // <<"set xlabel \"{/Symbol D}{/Symbol F}_{pj_2}\""<<endl
  	     // <<"plot 	\"deltaphi_pj2\" using 1:7"<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": p_{T,j_1}/p_{T,j_2}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/d(p_{T,j_1}:p_{T,j_2}) [fb]\""<<endl
  	     // <<"set xlabel \"p_{T,j_1}:p_{T,j_2}\""<<endl
  	     // <<"plot 	\"pt_j1:pt_j2\" using 1:7 title \"jet level\",\\"<<endl
  	     // <<"	\"pt_p1:pt_p2\" using 1:7 title \"associated partons\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": R_{pj}\""<<endl
  	     // <<"set ylabel \"d{/Symbol s}/dR_{pj} [fb]\""<<endl
  	     // <<"set xlabel \"R_{pj}\""<<endl
  	     // <<"plot 	\"R_pj1\" using 1:7 title \"hardest jet\",\\"<<endl
  	     // <<"	\"R_pj2\" using 1:7 title \"2nd hardest jet\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<": weights\""<<endl
  	     // <<"set ylabel \"Frequency\""<<endl
  	     // <<"set xlabel \"weight\""<<endl
  	     // <<"plot 	\"weights\" using 1:7 title \"weights\""<<endl
  	     // <<endl
  	     // <<"set style data fsteps"<<endl
  	     // <<"set title \""<<folder<<": Correlation plot\""<<endl
  	     // <<"set ylabel \"relative frequency\""<<endl
  	     // <<"set xlabel \"jet #\""<<endl
  	     // <<"set xrange [-0.5:8.5]"<<endl
  	     // <<"plot 	\"n_p1jnCorrelation\" using 2:4 title \"correlated jet of hardest parton\",\\"<<endl
  	     // <<"	\"n_p2jnCorrelation\" using 2:4 title \"correlated jet of 2nd hardest parton\""<<endl
  	     // <<endl
  	     // <<"set title \""<<folder<<" Correlation samples\""<<endl
  	     // <<"set xrange [0.0:5.0]"<<endl
  	     // <<"set ylabel \"# events\""<<endl
  	     // <<"set xtics (\"j1=p1 & j2=p2\" 0.5, \"j1=p2 & j2=p1\" 1.5, \"j1=p1 & j2!=p2\" 2.5, \"j1!=p1 & j2=p2\" 3.5, \"j1!=p1/2 & j2=p1/2\" 4.5)    "<<endl
  	     // <<"plot \"-\""<<endl
  	     // <<"0   0"<<endl
  	     // <<"1   "<<correlarray[0][0]<<endl
  	     // <<"2   "<<correlarray[0][2]<<endl
  	     // <<"3   "<<correlarray[0][1]<<endl
  	     // <<"4   "<<correlarray[1][0]<<endl
  	     // <<"5   "<<correlarray[1][1]<<endl
  	     // <<endl;

  if(strcmp(folder,"") != 0){
    if(chdir("..") != 0)
      cerr << "BasicHistogramCollegtion could not move into top directory.\n" << flush;
  }
}

/*
void BasicHistogramCollection::fill(const fastjet::PseudoJet & jet,const vector<fastjet::PseudoJet> & partons, int jetnumber){


  (*histo[jetnumber+1])+=jets[jetnumber].perp();
  (*histo[jetnumber+6])+=jets[jetnumber].eta();

  
  double Rn1=sqrt(pow(deltaphi(jets[jetnumber].phi_std(),partons[0].phi_std()),2)+pow(jets[jetnumber].eta()-partons[0].eta(),2));
  double Rn2=sqrt(pow(deltaphi(jets[jetnumber].phi_std(),partons[1].phi_std()),2)+pow(jets[jetnumber].eta()-partons[1].eta(),2));
  (*histo[jetnumber+8])+=min(R11,R12);

  (*histo[11])+=partons[0].perp();
  (*histo[12])+=partons[1].perp();
  (*histo[13])+=abs(partons[0].eta()-partons[1].eta());
  fastjet::PseudoJet sump = partons[0];
  sump+=partons[1];
  (*histo[14])+=sump.m();
  if(partons[0].eta()>partons[1].eta()){
    (*histo[15])+=deltaphi(partons[1].phi_std(),partons[0].phi_std());
  } 
  else{
    (*histo[15])+=deltaphi(partons[0].phi_std(),partons[1].phi_std());
  }
  (*histo[16])+=partons[0].eta();
  (*histo[17])+=partons[1].eta();
}
*/




SingleJetHistogramCollection::SingleJetHistogramCollection(){
  int n=100;

  histo[0] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  histo[1] =arHistogramPtr::Create(arHistogram(-7,7,n,true));
  histo[2] =arHistogramPtr::Create(arHistogram(0,10,n,true));
  histo[3] =arHistogramPtr::Create(arHistogram(0,300,n,true));
  histo[4] =arHistogramPtr::Create(arHistogram(-7,7,n,true));

}
SingleJetHistogramCollection::~SingleJetHistogramCollection(){
  // for (int i=0; i <= 4; i++){
  //   delete histo[i];
  // }
}


void SingleJetHistogramCollection::writetofolder(const char * folder){
  ofstream fileoutput[5];
  ofstream gnuplotfile;
  if(strcmp(folder,"") != 0){
    mkdir(folder,0777);
    if(chdir(folder) != 0)
      cerr << "BasicHistogramCollegtion could not move into subdirectory.\n" << flush;
  }

  fileoutput[0].open("pt_j");
  fileoutput[1].open("eta_j");
  fileoutput[2].open("R_pj");
  fileoutput[3].open("pt_p");
  fileoutput[4].open("eta_p");

  for (int i=0; i<=4; i++){
    (*histo[i]).simpleOutput(fileoutput[i],true,true,true);
    fileoutput[i].close();
  }

  gnuplotfile.open("SingleJetHistogramCollection.gp");
  gnuplotfile<<"set output \""<<folder<<".ps\""<<endl
	     <<"set terminal postscript enhanced color"<<endl
	     <<"#set style data boxes"<<endl
	     <<"#set style data histeps"<<endl
	     <<"set style data steps"<<endl
	     <<"set key right top"<<endl
	     <<endl
	     <<"set title \""<<folder<<": p_T of jet\""<<endl
	     <<"set ylabel \"d{/Symbol s}/dp_{T,j} [fb/GeV]\""<<endl
	     <<"set xlabel \"p_{T,j} [GeV]\""<<endl
	     <<"plot 	\"pt_j\" using 1:6 title \"jet level\",\\"<<endl
	     <<"	\"pt_p\" using 1:6 title \"associated partons\""<<endl
	     <<endl
	     <<"set title \""<<folder<<": rapidity of jet {/Symbol h}_{j}\""<<endl
	     <<"set ylabel \"d{/Symbol s}/d{/Symbol h}_{j} [fb]\""<<endl
	     <<"set xlabel \"{/Symbol h}_{j}\""<<endl
	     <<"plot 	\"eta_j\" using 1:6 title \"jet level\",\\"<<endl
	     <<"	\"eta_p\" using 1:6 title \"associated partons\""<<endl
	     <<endl
	     <<"set title \""<<folder<<": R_{pj}\""<<endl
	     <<"set ylabel \"d{/Symbol s}/dR_{pj} [fb]\""<<endl
	     <<"set xlabel \"R_{pj}\""<<endl
	     <<"plot 	\"R_pj\" using 1:6"<<endl;

  if(strcmp(folder,"") != 0){
    if(chdir("..") != 0)
      cerr << "SingleJetHistogramCollegtion could not move into top directory.\n" << flush;
  }
}

void SingleJetHistogramCollection::fill(const fastjet::PseudoJet & jet, const fastjet::PseudoJet & parton, double w){
  (*histo[0]).addWeighted(jet.perp(),w);
  (*histo[1]).addWeighted(jet.eta(),w);

  double Rpj=sqrt(pow(deltaphi(jet.phi_std(),parton.phi_std()),2)+pow(jet.eta()-parton.eta(),2));
  (*histo[2]).addWeighted(Rpj,w);

  (*histo[3]).addWeighted(parton.perp(),w);
  (*histo[4]).addWeighted(parton.eta(),w);

}


