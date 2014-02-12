#include "EventProxyOTF.h"


EventProxyOTF::EventProxyOTF(){}

EventProxyOTF::~EventProxyOTF(){}



void EventProxyOTF::init(std::vector<std::string> const& iFileNames){

	std::cout<<"EventProxyOTF "<<__func__<<std::endl;
	treeName_ = "efficiencyTree";

	EventProxyBase::init(iFileNames);

	fChain->SetBranchStatus("*",1);

	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Events_fUniqueID);
	fChain->SetBranchAddress("fBits", &fBits, &b_Events_fBits);
	fChain->SetBranchAddress("weight", &weight, &b_Events_weight);
	fChain->SetBranchAddress("pt", &pt, &b_Events_pt);
	fChain->SetBranchAddress("eta", &eta, &b_Events_eta);
	fChain->SetBranchAddress("phi", &phi, &b_Events_phi);
	fChain->SetBranchAddress("phiHit", &phiHit, &b_Events_phiHit);
	fChain->SetBranchAddress("etaHit", &etaHit, &b_Events_etaHit);
	fChain->SetBranchAddress("charge", &charge, &b_Events_charge);
	fChain->SetBranchAddress("l1ObjectsOtf", &l1ObjectsOtf_, &b_Events_l1ObjectsOtf_);
	fChain->SetBranchAddress("l1ObjectsOtf.fUniqueID", l1ObjectsOtf_fUniqueID, &b_l1ObjectsOtf_fUniqueID);
	fChain->SetBranchAddress("l1ObjectsOtf.fBits", l1ObjectsOtf_fBits, &b_l1ObjectsOtf_fBits);
	fChain->SetBranchAddress("l1ObjectsOtf.pt", l1ObjectsOtf_pt, &b_l1ObjectsOtf_pt);
	fChain->SetBranchAddress("l1ObjectsOtf.eta", l1ObjectsOtf_eta, &b_l1ObjectsOtf_eta);
	fChain->SetBranchAddress("l1ObjectsOtf.phi", l1ObjectsOtf_phi, &b_l1ObjectsOtf_phi);
	fChain->SetBranchAddress("l1ObjectsOtf.bx", l1ObjectsOtf_bx, &b_l1ObjectsOtf_bx);
	fChain->SetBranchAddress("l1ObjectsOtf.q", l1ObjectsOtf_q, &b_l1ObjectsOtf_q);
	fChain->SetBranchAddress("l1ObjectsOtf.charge", l1ObjectsOtf_charge, &b_l1ObjectsOtf_charge);
	fChain->SetBranchAddress("l1ObjectsOtf.type", l1ObjectsOtf_type, &b_l1ObjectsOtf_type);
	fChain->SetBranchAddress("l1ObjectsGmt", &l1ObjectsGmt_, &b_Events_l1ObjectsGmt_);
	fChain->SetBranchAddress("l1ObjectsGmt.fUniqueID", l1ObjectsGmt_fUniqueID, &b_l1ObjectsGmt_fUniqueID);
	fChain->SetBranchAddress("l1ObjectsGmt.fBits", l1ObjectsGmt_fBits, &b_l1ObjectsGmt_fBits);
	fChain->SetBranchAddress("l1ObjectsGmt.pt", l1ObjectsGmt_pt, &b_l1ObjectsGmt_pt);
	fChain->SetBranchAddress("l1ObjectsGmt.eta", l1ObjectsGmt_eta, &b_l1ObjectsGmt_eta);
	fChain->SetBranchAddress("l1ObjectsGmt.phi", l1ObjectsGmt_phi, &b_l1ObjectsGmt_phi);
	fChain->SetBranchAddress("l1ObjectsGmt.bx", l1ObjectsGmt_bx, &b_l1ObjectsGmt_bx);
	fChain->SetBranchAddress("l1ObjectsGmt.q", l1ObjectsGmt_q, &b_l1ObjectsGmt_q);
	fChain->SetBranchAddress("l1ObjectsGmt.charge", l1ObjectsGmt_charge, &b_l1ObjectsGmt_charge);
	fChain->SetBranchAddress("l1ObjectsGmt.type", l1ObjectsGmt_type, &b_l1ObjectsGmt_type);

}
