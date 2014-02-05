#ifndef RootAnalysis_ObjectMessenger_H
#define RootAnalysis_ObjectMessenger_H

#include <string>
#include <map>

class ObjectMessenger{

 public:

  ObjectMessenger(const std::string & aName);

  virtual ~ObjectMessenger(){;}
  
  const std::string & name(){return myName_;};

  ///Method for putting a object into messenger
  template<class T> void putObject(const T* obj, const std::string & flavour)= 0;

  ///Method for accessing object from the messenger
  template<class T> const T* getObject(const std::string & flavour) = 0;

  ///Method cleaning the messenger contents
  void clear();

 private:

  std::string myName_;

  ///Containers for carried objects to be defined in
  ///class

 
};

#endif

