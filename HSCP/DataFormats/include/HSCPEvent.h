#ifndef WarsawAnalysis_HSCPDataFormats_HSCPEvent_h
#define WarsawAnalysis_HSCPDataFormats_HSCPEvent_h

#include "TLorentzVector.h"
#include <map>
#include <vector>

#define MAX_CAND 16

class HSCPEvent;

class HSCPParticle{

  friend class HSCPEvent;
  friend class EventProxyHSCP;

  public:

  HSCPParticle(const HSCPEvent *aEvent);

  ~HSCPParticle(){}

  void clear();

  void init(int index);

  const float & getPt() const {return pt;}

  const float & getEta() const {return eta;}

  const float & getPhi() const {return phi;}

 protected:

  const HSCPEvent *myEvent;
  float pt, eta, phi;
  
};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HSCPEvent{

  friend class EventProxyHSCP;
  friend class HSCPParticle;

 public:

  HSCPEvent();

  ~HSCPEvent(){}

  ///Data member setters.
  void setRun(unsigned int x){runId = x;}

  void setEvent(unsigned long int x){eventId = x;}

  ///Reset class data members
  void clear();

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}

  int getNumberOfCandidates() const {return numberOfCandidates;}

  HSCPParticle getCandidate(unsigned int i) const;

 protected:

  ///Event run and number
  int runId;
  int eventId;
  int numberOfCandidates;

  Double_t pt[MAX_CAND], eta[MAX_CAND], phi[MAX_CAND]; 

};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
#endif
