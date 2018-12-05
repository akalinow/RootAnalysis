enum class PropertyEnum {
PDGId = 0, 
charge = 1, 
decayMode = 2, 
discriminator = 3, 
muonID = 4, 
typeOfMuon = 5, 
byCombinedIsolationDeltaBetaCorrRaw3Hits = 6, 
byIsolationMVArun2v1DBoldDMwLTraw = 7, 
byIsolationMVArun2v1DBoldDMwLTraw2017v2 = 8, 
byIsolationMVArun2v1DBnewDMwLTraw2017v2 = 9, 
dxy = 10, 
dz = 11, 
SIP = 12, 
tauID = 13, 
combreliso = 14, 
leadChargedParticlePt = 15, 
isGoodTriggerType = 16, 
FilterFired = 17, 
L3FilterFired = 18, 
L3FilterFiredLast = 19, 
mc_match = 20, 
rawPt = 21, 
area = 22, 
PUJetID = 23, 
jecUnc = 24, 
Flavour = 25, 
bDiscriminator = 26, 
bCSVscore = 27, 
PFjetID = 28, 
deepTau2017v1tauVSe = 29, 
deepTau2017v1tauVSmu = 30, 
deepTau2017v1tauVSjet = 31, 
deepTau2017v1tauVSall = 32, 
DPFTau_2016_v0tauVSall = 33, 
DPFTau_2016_v1tauVSall = 34, 
chargedIsoPtSum = 35, 
neutralIsoPtSum = 36, 
puCorrPtSum = 37, 
photonPtSumOutsideSignalCone = 38, 
nPhoton = 39, 
ptWeightedDetaStrip = 40, 
ptWeightedDphiStrip = 41, 
ptWeightedDrSignal = 42, 
ptWeightedDrIsolation = 43, 
eRatio = 44, 
dxy_Sig = 45, 
ip3d = 46, 
hasSecondaryVertex = 47, 
decayDistMag = 48, 
flightLengthSig = 49, 
gjAngleDiff = 50, 
NONE = 51
};

#include <map>
class PropertyEnumString { 
public:
static const std::map<std::string, PropertyEnum> enumMap;
};
