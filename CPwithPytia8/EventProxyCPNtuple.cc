#include <iostream>

#include "EventProxyCPNtuple.h"

EventProxyCPNtuple::EventProxyCPNtuple(){}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
EventProxyCPNtuple::~EventProxyCPNtuple(){}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
EventProxyBase* EventProxyCPNtuple::clone() const{

  return new EventProxyCPNtuple();
  
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void EventProxyCPNtuple::init(std::vector<std::string> const& iFileNames){

        treeName_ = "genAna/hTTCPTree";

	EventProxyBase::init(iFileNames);

	event = new HTTEvent();

	fChain->SetBranchAddress("p4Sum.",&event->genEvent_.p4SumPtr_, &b_p4Sum);
	
	fChain->SetBranchAddress("metNu.", &event->genEvent_.metNuPtr_, &b_metNu);
	fChain->SetBranchAddress("met.", &event->genEvent_.metPtr_, &b_met);
	fChain->SetBranchAddress("piMinus.", &event->genEvent_.piMinusPtr_, &b_piMinus);
	fChain->SetBranchAddress("piPlus.", &event->genEvent_.piPlusPtr_, &b_piPlus);
	fChain->SetBranchAddress("tauMinus.", &event->genEvent_.tauMinusPtr_, &b_tauMinus);
	fChain->SetBranchAddress("tauPlus.", &event->genEvent_.tauPlusPtr_, &b_tauPlus);
	fChain->SetBranchAddress("visTauMinus.", &event->genEvent_.visTauMinusPtr_, &b_visTauMinus);
	fChain->SetBranchAddress("visTauPlus.", &event->genEvent_.visTauPlusPtr_, &b_visTauPlus);
	fChain->SetBranchAddress("toyPiMinus.", &event->toyEvent_.piMinusPtr_, &b_toyPiMinus);
	fChain->SetBranchAddress("toyPiPlus.", &event->toyEvent_.piPlusPtr_, &b_toyPiPlus);
	fChain->SetBranchAddress("toyNeutralMinus.", &event->toyEvent_.neutralMinusPtr_, &b_toyNeutralMinus);
	fChain->SetBranchAddress("toyNeutralPlus.", &event->toyEvent_.neutralPlusPtr_, &b_toyNeutralPlus);
	fChain->SetBranchAddress("toyPiZeroMinus.", &event->toyEvent_.piZeroMinusPtr_, &b_toyPiZeroMinus);
	fChain->SetBranchAddress("toyPiZeroPlus.", &event->toyEvent_.piZeroPlusPtr_, &b_toyPiZeroPlus);
	fChain->SetBranchAddress("toyTauMinus.", &event->toyEvent_.tauMinusPtr_, &b_toyTauMinus);
	fChain->SetBranchAddress("toyTauPlus.", &event->toyEvent_.tauPlusPtr_, &b_toyTauPlus);
	//
	fChain->SetBranchAddress("bosonId", &event->bosonId_, &b_bosonId);
	fChain->SetBranchAddress("decModeMinus", &event->genEvent_.decModeMinus_, &b_decModeMinus);
	fChain->SetBranchAddress("decModePlus", &event->genEvent_.decModePlus_, &b_decModePlus);
	fChain->SetBranchAddress("toyDecModeMinus", &event->toyEvent_.decModeMinus_, &b_toyDecModeMinus);
	fChain->SetBranchAddress("toyDecModePlus", &event->toyEvent_.decModePlus_, &b_toyDecModePlus);
	fChain->SetBranchAddress("toyNChargedMinus", &event->toyEvent_.nChargedMinus_, &b_toyNChargedMinus);
	fChain->SetBranchAddress("toyNNeutralMinus", &event->toyEvent_.nNeutralMinus_, &b_toyNNeutralMinus);
	fChain->SetBranchAddress("toyNPiZeroMinus", &event->toyEvent_.nPiZeroMinus_, &b_toyNPiZeroMinus);
	fChain->SetBranchAddress("toyNChargedPlus", &event->toyEvent_.nChargedPlus_, &b_toyNChargedPlus);
	fChain->SetBranchAddress("toyNNeutralPlus", &event->toyEvent_.nNeutralPlus_, &b_toyNNeutralPlus);
	fChain->SetBranchAddress("toyNPiZeroPlus", &event->toyEvent_.nPiZeroPlus_, &b_toyNPiZeroPlus);
	//
	fChain->SetBranchAddress("thePV.", &event->genEvent_.thePVPtr_, &b_thePV);
	fChain->SetBranchAddress("svMinus.", &event->genEvent_.svMinusPtr_, &b_svMinus);
	fChain->SetBranchAddress("svPlus.", &event->genEvent_.svPlusPtr_, &b_svPlus);
	fChain->SetBranchAddress("nPiMinus.", &event->genEvent_.nPiMinusPtr_, &b_nPiMinus);
	fChain->SetBranchAddress("nPiPlus.", &event->genEvent_.nPiPlusPtr_, &b_nPiPlus);

	fChain->SetBranchAddress("toySvMinus.", &event->toyEvent_.svMinusPtr_, &b_toySvMinus);
	fChain->SetBranchAddress("toySvPlus.", &event->toyEvent_.svPlusPtr_, &b_toySvPlus);
	fChain->SetBranchAddress("toyNPiMinus.", &event->toyEvent_.nPiMinusPtr_, &b_toyNPiMinus);
	fChain->SetBranchAddress("toyNPiPlus.", &event->toyEvent_.nPiPlusPtr_, &b_toyNPiPlus);
	//
	fChain->SetBranchAddress("phi", &event->genEvent_.phiPtr_, &b_phi);
	fChain->SetBranchAddress("rho", &event->genEvent_.rhoPtr_, &b_rho);
	fChain->SetBranchAddress("phi2", &event->genEvent_.phi2Ptr_, &b_phi2);
	fChain->SetBranchAddress("phiRho", &event->genEvent_.phiRhoPtr_, &b_phiRho);
	//
	fChain->SetBranchAddress("yMinus", &event->genEvent_.yMinus_, &b_yMinus);
	fChain->SetBranchAddress("yPlus", &event->genEvent_.yPlus_, &b_yPlus);
	fChain->SetBranchAddress("yMinus2", &event->genEvent_.yMinus2_, &b_yMinus2);
	fChain->SetBranchAddress("yPlus2", &event->genEvent_.yPlus2_, &b_yPlus2);
	fChain->SetBranchAddress("yMinusLab", &event->genEvent_.yMinusLab_, &b_yMinusLab);
	fChain->SetBranchAddress("yPlusLab", &event->genEvent_.yPlusLab_, &b_yPlusLab);
	fChain->SetBranchAddress("yMinusLab2", &event->genEvent_.yMinusLab2_, &b_yMinusLab2);
	fChain->SetBranchAddress("yPlusLab2", &event->genEvent_.yPlusLab2_, &b_yPlusLab2);
	fChain->SetBranchAddress("yToyMinus", &event->toyEvent_.yMinus_, &b_yToyMinus);
	fChain->SetBranchAddress("yToyPlus", &event->toyEvent_.yPlus_, &b_yToyPlus);
	fChain->SetBranchAddress("isoMinus", &event->toyEvent_.isoMinus_, &b_isoMinus);
	fChain->SetBranchAddress("isoPlus", &event->toyEvent_.isoPlus_, &b_isoPlus);
	fChain->SetBranchAddress("outerMinus", &event->toyEvent_.outerMinus_, &b_outerMinus);
	fChain->SetBranchAddress("outerPlus", &event->toyEvent_.outerPlus_, &b_outerPlus);
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
