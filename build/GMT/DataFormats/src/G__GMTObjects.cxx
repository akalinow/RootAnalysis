// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__GMTObjects
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "EventObj.h"
#include "GenObj.h"
#include "GenObjColl.h"
#include "L1Obj.h"
#include "L1ObjColl.h"
#include "L1PhaseIIObj.h"
#include "L1PhaseIIObjColl.h"
#include "RecoMuon.h"
#include "RecoMuonObj.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_EventObj(void *p = 0);
   static void *newArray_EventObj(Long_t size, void *p);
   static void delete_EventObj(void *p);
   static void deleteArray_EventObj(void *p);
   static void destruct_EventObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventObj*)
   {
      ::EventObj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventObj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EventObj", ::EventObj::Class_Version(), "EventObj.h", 8,
                  typeid(::EventObj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventObj::Dictionary, isa_proxy, 4,
                  sizeof(::EventObj) );
      instance.SetNew(&new_EventObj);
      instance.SetNewArray(&newArray_EventObj);
      instance.SetDelete(&delete_EventObj);
      instance.SetDeleteArray(&deleteArray_EventObj);
      instance.SetDestructor(&destruct_EventObj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventObj*)
   {
      return GenerateInitInstanceLocal((::EventObj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EventObj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GenObj(void *p = 0);
   static void *newArray_GenObj(Long_t size, void *p);
   static void delete_GenObj(void *p);
   static void deleteArray_GenObj(void *p);
   static void destruct_GenObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GenObj*)
   {
      ::GenObj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GenObj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GenObj", ::GenObj::Class_Version(), "GenObj.h", 8,
                  typeid(::GenObj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GenObj::Dictionary, isa_proxy, 4,
                  sizeof(::GenObj) );
      instance.SetNew(&new_GenObj);
      instance.SetNewArray(&newArray_GenObj);
      instance.SetDelete(&delete_GenObj);
      instance.SetDeleteArray(&deleteArray_GenObj);
      instance.SetDestructor(&destruct_GenObj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GenObj*)
   {
      return GenerateInitInstanceLocal((::GenObj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GenObj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GenObjColl(void *p = 0);
   static void *newArray_GenObjColl(Long_t size, void *p);
   static void delete_GenObjColl(void *p);
   static void deleteArray_GenObjColl(void *p);
   static void destruct_GenObjColl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GenObjColl*)
   {
      ::GenObjColl *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GenObjColl >(0);
      static ::ROOT::TGenericClassInfo 
         instance("GenObjColl", ::GenObjColl::Class_Version(), "GenObjColl.h", 10,
                  typeid(::GenObjColl), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::GenObjColl::Dictionary, isa_proxy, 4,
                  sizeof(::GenObjColl) );
      instance.SetNew(&new_GenObjColl);
      instance.SetNewArray(&newArray_GenObjColl);
      instance.SetDelete(&delete_GenObjColl);
      instance.SetDeleteArray(&deleteArray_GenObjColl);
      instance.SetDestructor(&destruct_GenObjColl);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::GenObjColl*)
   {
      return GenerateInitInstanceLocal((::GenObjColl*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GenObjColl*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_L1Obj(void *p = 0);
   static void *newArray_L1Obj(Long_t size, void *p);
   static void delete_L1Obj(void *p);
   static void deleteArray_L1Obj(void *p);
   static void destruct_L1Obj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1Obj*)
   {
      ::L1Obj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::L1Obj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("L1Obj", ::L1Obj::Class_Version(), "L1Obj.h", 15,
                  typeid(::L1Obj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::L1Obj::Dictionary, isa_proxy, 4,
                  sizeof(::L1Obj) );
      instance.SetNew(&new_L1Obj);
      instance.SetNewArray(&newArray_L1Obj);
      instance.SetDelete(&delete_L1Obj);
      instance.SetDeleteArray(&deleteArray_L1Obj);
      instance.SetDestructor(&destruct_L1Obj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::L1Obj*)
   {
      return GenerateInitInstanceLocal((::L1Obj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1Obj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_L1ObjColl(void *p = 0);
   static void *newArray_L1ObjColl(Long_t size, void *p);
   static void delete_L1ObjColl(void *p);
   static void deleteArray_L1ObjColl(void *p);
   static void destruct_L1ObjColl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1ObjColl*)
   {
      ::L1ObjColl *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::L1ObjColl >(0);
      static ::ROOT::TGenericClassInfo 
         instance("L1ObjColl", ::L1ObjColl::Class_Version(), "L1ObjColl.h", 9,
                  typeid(::L1ObjColl), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::L1ObjColl::Dictionary, isa_proxy, 4,
                  sizeof(::L1ObjColl) );
      instance.SetNew(&new_L1ObjColl);
      instance.SetNewArray(&newArray_L1ObjColl);
      instance.SetDelete(&delete_L1ObjColl);
      instance.SetDeleteArray(&deleteArray_L1ObjColl);
      instance.SetDestructor(&destruct_L1ObjColl);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::L1ObjColl*)
   {
      return GenerateInitInstanceLocal((::L1ObjColl*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1ObjColl*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_L1PhaseIIObj(void *p = 0);
   static void *newArray_L1PhaseIIObj(Long_t size, void *p);
   static void delete_L1PhaseIIObj(void *p);
   static void deleteArray_L1PhaseIIObj(void *p);
   static void destruct_L1PhaseIIObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1PhaseIIObj*)
   {
      ::L1PhaseIIObj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::L1PhaseIIObj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("L1PhaseIIObj", ::L1PhaseIIObj::Class_Version(), "L1PhaseIIObj.h", 8,
                  typeid(::L1PhaseIIObj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::L1PhaseIIObj::Dictionary, isa_proxy, 4,
                  sizeof(::L1PhaseIIObj) );
      instance.SetNew(&new_L1PhaseIIObj);
      instance.SetNewArray(&newArray_L1PhaseIIObj);
      instance.SetDelete(&delete_L1PhaseIIObj);
      instance.SetDeleteArray(&deleteArray_L1PhaseIIObj);
      instance.SetDestructor(&destruct_L1PhaseIIObj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::L1PhaseIIObj*)
   {
      return GenerateInitInstanceLocal((::L1PhaseIIObj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1PhaseIIObj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *L1PhaseIIObjColl_Dictionary();
   static void L1PhaseIIObjColl_TClassManip(TClass*);
   static void *new_L1PhaseIIObjColl(void *p = 0);
   static void *newArray_L1PhaseIIObjColl(Long_t size, void *p);
   static void delete_L1PhaseIIObjColl(void *p);
   static void deleteArray_L1PhaseIIObjColl(void *p);
   static void destruct_L1PhaseIIObjColl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1PhaseIIObjColl*)
   {
      ::L1PhaseIIObjColl *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::L1PhaseIIObjColl));
      static ::ROOT::TGenericClassInfo 
         instance("L1PhaseIIObjColl", "L1PhaseIIObjColl.h", 8,
                  typeid(::L1PhaseIIObjColl), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &L1PhaseIIObjColl_Dictionary, isa_proxy, 4,
                  sizeof(::L1PhaseIIObjColl) );
      instance.SetNew(&new_L1PhaseIIObjColl);
      instance.SetNewArray(&newArray_L1PhaseIIObjColl);
      instance.SetDelete(&delete_L1PhaseIIObjColl);
      instance.SetDeleteArray(&deleteArray_L1PhaseIIObjColl);
      instance.SetDestructor(&destruct_L1PhaseIIObjColl);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::L1PhaseIIObjColl*)
   {
      return GenerateInitInstanceLocal((::L1PhaseIIObjColl*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1PhaseIIObjColl*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *L1PhaseIIObjColl_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::L1PhaseIIObjColl*)0x0)->GetClass();
      L1PhaseIIObjColl_TClassManip(theClass);
   return theClass;
   }

   static void L1PhaseIIObjColl_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_RecoMuonObj(void *p = 0);
   static void *newArray_RecoMuonObj(Long_t size, void *p);
   static void delete_RecoMuonObj(void *p);
   static void deleteArray_RecoMuonObj(void *p);
   static void destruct_RecoMuonObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RecoMuonObj*)
   {
      ::RecoMuonObj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RecoMuonObj >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RecoMuonObj", ::RecoMuonObj::Class_Version(), "RecoMuonObj.h", 7,
                  typeid(::RecoMuonObj), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RecoMuonObj::Dictionary, isa_proxy, 4,
                  sizeof(::RecoMuonObj) );
      instance.SetNew(&new_RecoMuonObj);
      instance.SetNewArray(&newArray_RecoMuonObj);
      instance.SetDelete(&delete_RecoMuonObj);
      instance.SetDeleteArray(&deleteArray_RecoMuonObj);
      instance.SetDestructor(&destruct_RecoMuonObj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RecoMuonObj*)
   {
      return GenerateInitInstanceLocal((::RecoMuonObj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RecoMuonObj*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RecoMuon(void *p = 0);
   static void *newArray_RecoMuon(Long_t size, void *p);
   static void delete_RecoMuon(void *p);
   static void deleteArray_RecoMuon(void *p);
   static void destruct_RecoMuon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RecoMuon*)
   {
      ::RecoMuon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RecoMuon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RecoMuon", ::RecoMuon::Class_Version(), "RecoMuon.h", 10,
                  typeid(::RecoMuon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RecoMuon::Dictionary, isa_proxy, 4,
                  sizeof(::RecoMuon) );
      instance.SetNew(&new_RecoMuon);
      instance.SetNewArray(&newArray_RecoMuon);
      instance.SetDelete(&delete_RecoMuon);
      instance.SetDeleteArray(&deleteArray_RecoMuon);
      instance.SetDestructor(&destruct_RecoMuon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RecoMuon*)
   {
      return GenerateInitInstanceLocal((::RecoMuon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RecoMuon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EventObj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventObj::Class_Name()
{
   return "EventObj";
}

//______________________________________________________________________________
const char *EventObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GenObj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenObj::Class_Name()
{
   return "GenObj";
}

//______________________________________________________________________________
const char *GenObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GenObjColl::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenObjColl::Class_Name()
{
   return "GenObjColl";
}

//______________________________________________________________________________
const char *GenObjColl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenObjColl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenObjColl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenObjColl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr L1Obj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *L1Obj::Class_Name()
{
   return "L1Obj";
}

//______________________________________________________________________________
const char *L1Obj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int L1Obj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1Obj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1Obj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr L1ObjColl::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *L1ObjColl::Class_Name()
{
   return "L1ObjColl";
}

//______________________________________________________________________________
const char *L1ObjColl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int L1ObjColl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1ObjColl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1ObjColl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr L1PhaseIIObj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *L1PhaseIIObj::Class_Name()
{
   return "L1PhaseIIObj";
}

//______________________________________________________________________________
const char *L1PhaseIIObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1PhaseIIObj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int L1PhaseIIObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1PhaseIIObj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1PhaseIIObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1PhaseIIObj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1PhaseIIObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1PhaseIIObj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RecoMuonObj::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RecoMuonObj::Class_Name()
{
   return "RecoMuonObj";
}

//______________________________________________________________________________
const char *RecoMuonObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoMuonObj*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RecoMuonObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoMuonObj*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RecoMuonObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoMuonObj*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RecoMuonObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoMuonObj*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RecoMuon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RecoMuon::Class_Name()
{
   return "RecoMuon";
}

//______________________________________________________________________________
const char *RecoMuon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoMuon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RecoMuon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoMuon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RecoMuon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoMuon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RecoMuon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoMuon*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EventObj::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventObj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(EventObj::Class(),this);
   } else {
      R__b.WriteClassBuffer(EventObj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventObj(void *p) {
      return  p ? new(p) ::EventObj : new ::EventObj;
   }
   static void *newArray_EventObj(Long_t nElements, void *p) {
      return p ? new(p) ::EventObj[nElements] : new ::EventObj[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventObj(void *p) {
      delete ((::EventObj*)p);
   }
   static void deleteArray_EventObj(void *p) {
      delete [] ((::EventObj*)p);
   }
   static void destruct_EventObj(void *p) {
      typedef ::EventObj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::EventObj

//______________________________________________________________________________
void GenObj::Streamer(TBuffer &R__b)
{
   // Stream an object of class GenObj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GenObj::Class(),this);
   } else {
      R__b.WriteClassBuffer(GenObj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GenObj(void *p) {
      return  p ? new(p) ::GenObj : new ::GenObj;
   }
   static void *newArray_GenObj(Long_t nElements, void *p) {
      return p ? new(p) ::GenObj[nElements] : new ::GenObj[nElements];
   }
   // Wrapper around operator delete
   static void delete_GenObj(void *p) {
      delete ((::GenObj*)p);
   }
   static void deleteArray_GenObj(void *p) {
      delete [] ((::GenObj*)p);
   }
   static void destruct_GenObj(void *p) {
      typedef ::GenObj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GenObj

//______________________________________________________________________________
void GenObjColl::Streamer(TBuffer &R__b)
{
   // Stream an object of class GenObjColl.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(GenObjColl::Class(),this);
   } else {
      R__b.WriteClassBuffer(GenObjColl::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_GenObjColl(void *p) {
      return  p ? new(p) ::GenObjColl : new ::GenObjColl;
   }
   static void *newArray_GenObjColl(Long_t nElements, void *p) {
      return p ? new(p) ::GenObjColl[nElements] : new ::GenObjColl[nElements];
   }
   // Wrapper around operator delete
   static void delete_GenObjColl(void *p) {
      delete ((::GenObjColl*)p);
   }
   static void deleteArray_GenObjColl(void *p) {
      delete [] ((::GenObjColl*)p);
   }
   static void destruct_GenObjColl(void *p) {
      typedef ::GenObjColl current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::GenObjColl

//______________________________________________________________________________
void L1Obj::Streamer(TBuffer &R__b)
{
   // Stream an object of class L1Obj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(L1Obj::Class(),this);
   } else {
      R__b.WriteClassBuffer(L1Obj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_L1Obj(void *p) {
      return  p ? new(p) ::L1Obj : new ::L1Obj;
   }
   static void *newArray_L1Obj(Long_t nElements, void *p) {
      return p ? new(p) ::L1Obj[nElements] : new ::L1Obj[nElements];
   }
   // Wrapper around operator delete
   static void delete_L1Obj(void *p) {
      delete ((::L1Obj*)p);
   }
   static void deleteArray_L1Obj(void *p) {
      delete [] ((::L1Obj*)p);
   }
   static void destruct_L1Obj(void *p) {
      typedef ::L1Obj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::L1Obj

//______________________________________________________________________________
void L1ObjColl::Streamer(TBuffer &R__b)
{
   // Stream an object of class L1ObjColl.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(L1ObjColl::Class(),this);
   } else {
      R__b.WriteClassBuffer(L1ObjColl::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_L1ObjColl(void *p) {
      return  p ? new(p) ::L1ObjColl : new ::L1ObjColl;
   }
   static void *newArray_L1ObjColl(Long_t nElements, void *p) {
      return p ? new(p) ::L1ObjColl[nElements] : new ::L1ObjColl[nElements];
   }
   // Wrapper around operator delete
   static void delete_L1ObjColl(void *p) {
      delete ((::L1ObjColl*)p);
   }
   static void deleteArray_L1ObjColl(void *p) {
      delete [] ((::L1ObjColl*)p);
   }
   static void destruct_L1ObjColl(void *p) {
      typedef ::L1ObjColl current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::L1ObjColl

//______________________________________________________________________________
void L1PhaseIIObj::Streamer(TBuffer &R__b)
{
   // Stream an object of class L1PhaseIIObj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(L1PhaseIIObj::Class(),this);
   } else {
      R__b.WriteClassBuffer(L1PhaseIIObj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_L1PhaseIIObj(void *p) {
      return  p ? new(p) ::L1PhaseIIObj : new ::L1PhaseIIObj;
   }
   static void *newArray_L1PhaseIIObj(Long_t nElements, void *p) {
      return p ? new(p) ::L1PhaseIIObj[nElements] : new ::L1PhaseIIObj[nElements];
   }
   // Wrapper around operator delete
   static void delete_L1PhaseIIObj(void *p) {
      delete ((::L1PhaseIIObj*)p);
   }
   static void deleteArray_L1PhaseIIObj(void *p) {
      delete [] ((::L1PhaseIIObj*)p);
   }
   static void destruct_L1PhaseIIObj(void *p) {
      typedef ::L1PhaseIIObj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::L1PhaseIIObj

namespace ROOT {
   // Wrappers around operator new
   static void *new_L1PhaseIIObjColl(void *p) {
      return  p ? new(p) ::L1PhaseIIObjColl : new ::L1PhaseIIObjColl;
   }
   static void *newArray_L1PhaseIIObjColl(Long_t nElements, void *p) {
      return p ? new(p) ::L1PhaseIIObjColl[nElements] : new ::L1PhaseIIObjColl[nElements];
   }
   // Wrapper around operator delete
   static void delete_L1PhaseIIObjColl(void *p) {
      delete ((::L1PhaseIIObjColl*)p);
   }
   static void deleteArray_L1PhaseIIObjColl(void *p) {
      delete [] ((::L1PhaseIIObjColl*)p);
   }
   static void destruct_L1PhaseIIObjColl(void *p) {
      typedef ::L1PhaseIIObjColl current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::L1PhaseIIObjColl

//______________________________________________________________________________
void RecoMuonObj::Streamer(TBuffer &R__b)
{
   // Stream an object of class RecoMuonObj.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RecoMuonObj::Class(),this);
   } else {
      R__b.WriteClassBuffer(RecoMuonObj::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RecoMuonObj(void *p) {
      return  p ? new(p) ::RecoMuonObj : new ::RecoMuonObj;
   }
   static void *newArray_RecoMuonObj(Long_t nElements, void *p) {
      return p ? new(p) ::RecoMuonObj[nElements] : new ::RecoMuonObj[nElements];
   }
   // Wrapper around operator delete
   static void delete_RecoMuonObj(void *p) {
      delete ((::RecoMuonObj*)p);
   }
   static void deleteArray_RecoMuonObj(void *p) {
      delete [] ((::RecoMuonObj*)p);
   }
   static void destruct_RecoMuonObj(void *p) {
      typedef ::RecoMuonObj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RecoMuonObj

//______________________________________________________________________________
void RecoMuon::Streamer(TBuffer &R__b)
{
   // Stream an object of class RecoMuon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RecoMuon::Class(),this);
   } else {
      R__b.WriteClassBuffer(RecoMuon::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RecoMuon(void *p) {
      return  p ? new(p) ::RecoMuon : new ::RecoMuon;
   }
   static void *newArray_RecoMuon(Long_t nElements, void *p) {
      return p ? new(p) ::RecoMuon[nElements] : new ::RecoMuon[nElements];
   }
   // Wrapper around operator delete
   static void delete_RecoMuon(void *p) {
      delete ((::RecoMuon*)p);
   }
   static void deleteArray_RecoMuon(void *p) {
      delete [] ((::RecoMuon*)p);
   }
   static void destruct_RecoMuon(void *p) {
      typedef ::RecoMuon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RecoMuon

namespace ROOT {
   static TClass *vectorlElongsPdoublegR_Dictionary();
   static void vectorlElongsPdoublegR_TClassManip(TClass*);
   static void *new_vectorlElongsPdoublegR(void *p = 0);
   static void *newArray_vectorlElongsPdoublegR(Long_t size, void *p);
   static void delete_vectorlElongsPdoublegR(void *p);
   static void deleteArray_vectorlElongsPdoublegR(void *p);
   static void destruct_vectorlElongsPdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<long double>*)
   {
      vector<long double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<long double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<long double>", -2, "vector", 210,
                  typeid(vector<long double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElongsPdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<long double>) );
      instance.SetNew(&new_vectorlElongsPdoublegR);
      instance.SetNewArray(&newArray_vectorlElongsPdoublegR);
      instance.SetDelete(&delete_vectorlElongsPdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlElongsPdoublegR);
      instance.SetDestructor(&destruct_vectorlElongsPdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<long double> >()));

      ::ROOT::AddClassAlternate("vector<long double>","std::vector<long double, std::allocator<long double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<long double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElongsPdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<long double>*)0x0)->GetClass();
      vectorlElongsPdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElongsPdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElongsPdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long double> : new vector<long double>;
   }
   static void *newArray_vectorlElongsPdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long double>[nElements] : new vector<long double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElongsPdoublegR(void *p) {
      delete ((vector<long double>*)p);
   }
   static void deleteArray_vectorlElongsPdoublegR(void *p) {
      delete [] ((vector<long double>*)p);
   }
   static void destruct_vectorlElongsPdoublegR(void *p) {
      typedef vector<long double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<long double>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = 0);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 518,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));

      ::ROOT::AddClassAlternate("vector<bool>","std::vector<bool, std::allocator<bool> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<bool>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bool>*)0x0)->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete ((vector<bool>*)p);
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] ((vector<bool>*)p);
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace ROOT {
   static TClass *vectorlERecoMuonObjgR_Dictionary();
   static void vectorlERecoMuonObjgR_TClassManip(TClass*);
   static void *new_vectorlERecoMuonObjgR(void *p = 0);
   static void *newArray_vectorlERecoMuonObjgR(Long_t size, void *p);
   static void delete_vectorlERecoMuonObjgR(void *p);
   static void deleteArray_vectorlERecoMuonObjgR(void *p);
   static void destruct_vectorlERecoMuonObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RecoMuonObj>*)
   {
      vector<RecoMuonObj> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RecoMuonObj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RecoMuonObj>", -2, "vector", 210,
                  typeid(vector<RecoMuonObj>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERecoMuonObjgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RecoMuonObj>) );
      instance.SetNew(&new_vectorlERecoMuonObjgR);
      instance.SetNewArray(&newArray_vectorlERecoMuonObjgR);
      instance.SetDelete(&delete_vectorlERecoMuonObjgR);
      instance.SetDeleteArray(&deleteArray_vectorlERecoMuonObjgR);
      instance.SetDestructor(&destruct_vectorlERecoMuonObjgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RecoMuonObj> >()));

      ::ROOT::AddClassAlternate("vector<RecoMuonObj>","std::vector<RecoMuonObj, std::allocator<RecoMuonObj> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RecoMuonObj>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERecoMuonObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RecoMuonObj>*)0x0)->GetClass();
      vectorlERecoMuonObjgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERecoMuonObjgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERecoMuonObjgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RecoMuonObj> : new vector<RecoMuonObj>;
   }
   static void *newArray_vectorlERecoMuonObjgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RecoMuonObj>[nElements] : new vector<RecoMuonObj>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERecoMuonObjgR(void *p) {
      delete ((vector<RecoMuonObj>*)p);
   }
   static void deleteArray_vectorlERecoMuonObjgR(void *p) {
      delete [] ((vector<RecoMuonObj>*)p);
   }
   static void destruct_vectorlERecoMuonObjgR(void *p) {
      typedef vector<RecoMuonObj> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RecoMuonObj>

namespace ROOT {
   static TClass *vectorlEL1PhaseIIObjgR_Dictionary();
   static void vectorlEL1PhaseIIObjgR_TClassManip(TClass*);
   static void *new_vectorlEL1PhaseIIObjgR(void *p = 0);
   static void *newArray_vectorlEL1PhaseIIObjgR(Long_t size, void *p);
   static void delete_vectorlEL1PhaseIIObjgR(void *p);
   static void deleteArray_vectorlEL1PhaseIIObjgR(void *p);
   static void destruct_vectorlEL1PhaseIIObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<L1PhaseIIObj>*)
   {
      vector<L1PhaseIIObj> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<L1PhaseIIObj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<L1PhaseIIObj>", -2, "vector", 210,
                  typeid(vector<L1PhaseIIObj>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEL1PhaseIIObjgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<L1PhaseIIObj>) );
      instance.SetNew(&new_vectorlEL1PhaseIIObjgR);
      instance.SetNewArray(&newArray_vectorlEL1PhaseIIObjgR);
      instance.SetDelete(&delete_vectorlEL1PhaseIIObjgR);
      instance.SetDeleteArray(&deleteArray_vectorlEL1PhaseIIObjgR);
      instance.SetDestructor(&destruct_vectorlEL1PhaseIIObjgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<L1PhaseIIObj> >()));

      ::ROOT::AddClassAlternate("vector<L1PhaseIIObj>","std::vector<L1PhaseIIObj, std::allocator<L1PhaseIIObj> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<L1PhaseIIObj>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEL1PhaseIIObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<L1PhaseIIObj>*)0x0)->GetClass();
      vectorlEL1PhaseIIObjgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEL1PhaseIIObjgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEL1PhaseIIObjgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<L1PhaseIIObj> : new vector<L1PhaseIIObj>;
   }
   static void *newArray_vectorlEL1PhaseIIObjgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<L1PhaseIIObj>[nElements] : new vector<L1PhaseIIObj>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEL1PhaseIIObjgR(void *p) {
      delete ((vector<L1PhaseIIObj>*)p);
   }
   static void deleteArray_vectorlEL1PhaseIIObjgR(void *p) {
      delete [] ((vector<L1PhaseIIObj>*)p);
   }
   static void destruct_vectorlEL1PhaseIIObjgR(void *p) {
      typedef vector<L1PhaseIIObj> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<L1PhaseIIObj>

namespace ROOT {
   static TClass *vectorlEL1ObjgR_Dictionary();
   static void vectorlEL1ObjgR_TClassManip(TClass*);
   static void *new_vectorlEL1ObjgR(void *p = 0);
   static void *newArray_vectorlEL1ObjgR(Long_t size, void *p);
   static void delete_vectorlEL1ObjgR(void *p);
   static void deleteArray_vectorlEL1ObjgR(void *p);
   static void destruct_vectorlEL1ObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<L1Obj>*)
   {
      vector<L1Obj> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<L1Obj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<L1Obj>", -2, "vector", 210,
                  typeid(vector<L1Obj>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEL1ObjgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<L1Obj>) );
      instance.SetNew(&new_vectorlEL1ObjgR);
      instance.SetNewArray(&newArray_vectorlEL1ObjgR);
      instance.SetDelete(&delete_vectorlEL1ObjgR);
      instance.SetDeleteArray(&deleteArray_vectorlEL1ObjgR);
      instance.SetDestructor(&destruct_vectorlEL1ObjgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<L1Obj> >()));

      ::ROOT::AddClassAlternate("vector<L1Obj>","std::vector<L1Obj, std::allocator<L1Obj> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<L1Obj>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEL1ObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<L1Obj>*)0x0)->GetClass();
      vectorlEL1ObjgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEL1ObjgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEL1ObjgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<L1Obj> : new vector<L1Obj>;
   }
   static void *newArray_vectorlEL1ObjgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<L1Obj>[nElements] : new vector<L1Obj>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEL1ObjgR(void *p) {
      delete ((vector<L1Obj>*)p);
   }
   static void deleteArray_vectorlEL1ObjgR(void *p) {
      delete [] ((vector<L1Obj>*)p);
   }
   static void destruct_vectorlEL1ObjgR(void *p) {
      typedef vector<L1Obj> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<L1Obj>

namespace ROOT {
   static TClass *vectorlEGenObjgR_Dictionary();
   static void vectorlEGenObjgR_TClassManip(TClass*);
   static void *new_vectorlEGenObjgR(void *p = 0);
   static void *newArray_vectorlEGenObjgR(Long_t size, void *p);
   static void delete_vectorlEGenObjgR(void *p);
   static void deleteArray_vectorlEGenObjgR(void *p);
   static void destruct_vectorlEGenObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<GenObj>*)
   {
      vector<GenObj> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<GenObj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<GenObj>", -2, "vector", 210,
                  typeid(vector<GenObj>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEGenObjgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<GenObj>) );
      instance.SetNew(&new_vectorlEGenObjgR);
      instance.SetNewArray(&newArray_vectorlEGenObjgR);
      instance.SetDelete(&delete_vectorlEGenObjgR);
      instance.SetDeleteArray(&deleteArray_vectorlEGenObjgR);
      instance.SetDestructor(&destruct_vectorlEGenObjgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<GenObj> >()));

      ::ROOT::AddClassAlternate("vector<GenObj>","std::vector<GenObj, std::allocator<GenObj> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<GenObj>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEGenObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<GenObj>*)0x0)->GetClass();
      vectorlEGenObjgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEGenObjgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEGenObjgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GenObj> : new vector<GenObj>;
   }
   static void *newArray_vectorlEGenObjgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<GenObj>[nElements] : new vector<GenObj>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEGenObjgR(void *p) {
      delete ((vector<GenObj>*)p);
   }
   static void deleteArray_vectorlEGenObjgR(void *p) {
      delete [] ((vector<GenObj>*)p);
   }
   static void destruct_vectorlEGenObjgR(void *p) {
      typedef vector<GenObj> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<GenObj>

namespace {
  void TriggerDictionaryInitialization_libGMTObjects_Impl() {
    static const char* headers[] = {
"EventObj.h",
"GenObj.h",
"GenObjColl.h",
"L1Obj.h",
"L1ObjColl.h",
"L1PhaseIIObj.h",
"L1PhaseIIObjColl.h",
"RecoMuon.h",
"RecoMuonObj.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/eos/user/a/almuhamm/02.TriggerEff/RootAnalysis/GMT/DataFormats/src/../include",
"/eos/user/a/almuhamm/02.TriggerEff/RootAnalysis/GMT/DataFormats/src",
"/usr/include/root",
"/eos/home-a/almuhamm/02.TriggerEff/RootAnalysis/build/GMT/DataFormats/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libGMTObjects dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate("$clingAutoload$EventObj.h")))  EventObj;
class __attribute__((annotate("$clingAutoload$GenObj.h")))  GenObj;
class __attribute__((annotate("$clingAutoload$GenObjColl.h")))  GenObjColl;
class __attribute__((annotate("$clingAutoload$L1Obj.h")))  L1Obj;
class __attribute__((annotate("$clingAutoload$L1ObjColl.h")))  L1ObjColl;
class __attribute__((annotate("$clingAutoload$L1PhaseIIObj.h")))  L1PhaseIIObj;
class __attribute__((annotate("$clingAutoload$L1PhaseIIObjColl.h")))  L1PhaseIIObjColl;
class __attribute__((annotate("$clingAutoload$RecoMuonObj.h")))  __attribute__((annotate("$clingAutoload$RecoMuon.h")))  RecoMuonObj;
class __attribute__((annotate("$clingAutoload$RecoMuon.h")))  RecoMuon;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libGMTObjects dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "EventObj.h"
#include "GenObj.h"
#include "GenObjColl.h"
#include "L1Obj.h"
#include "L1ObjColl.h"
#include "L1PhaseIIObj.h"
#include "L1PhaseIIObjColl.h"
#include "RecoMuon.h"
#include "RecoMuonObj.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EventObj", payloadCode, "@",
"GenObj", payloadCode, "@",
"GenObjColl", payloadCode, "@",
"L1Obj", payloadCode, "@",
"L1ObjColl", payloadCode, "@",
"L1PhaseIIObj", payloadCode, "@",
"L1PhaseIIObjColl", payloadCode, "@",
"RecoMuon", payloadCode, "@",
"RecoMuonObj", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGMTObjects",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGMTObjects_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGMTObjects_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGMTObjects() {
  TriggerDictionaryInitialization_libGMTObjects_Impl();
}
