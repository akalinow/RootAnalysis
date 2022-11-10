// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__OMTFObjects
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

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_EventObj(void *p = nullptr);
   static void *newArray_EventObj(Long_t size, void *p);
   static void delete_EventObj(void *p);
   static void deleteArray_EventObj(void *p);
   static void destruct_EventObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventObj*)
   {
      ::EventObj *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventObj >(nullptr);
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
      return GenerateInitInstanceLocal((::EventObj*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EventObj*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GenObj(void *p = nullptr);
   static void *newArray_GenObj(Long_t size, void *p);
   static void delete_GenObj(void *p);
   static void deleteArray_GenObj(void *p);
   static void destruct_GenObj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GenObj*)
   {
      ::GenObj *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GenObj >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("GenObj", ::GenObj::Class_Version(), "GenObj.h", 6,
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
      return GenerateInitInstanceLocal((::GenObj*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GenObj*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_GenObjColl(void *p = nullptr);
   static void *newArray_GenObjColl(Long_t size, void *p);
   static void delete_GenObjColl(void *p);
   static void deleteArray_GenObjColl(void *p);
   static void destruct_GenObjColl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::GenObjColl*)
   {
      ::GenObjColl *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::GenObjColl >(nullptr);
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
      return GenerateInitInstanceLocal((::GenObjColl*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::GenObjColl*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_L1Obj(void *p = nullptr);
   static void *newArray_L1Obj(Long_t size, void *p);
   static void delete_L1Obj(void *p);
   static void deleteArray_L1Obj(void *p);
   static void destruct_L1Obj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1Obj*)
   {
      ::L1Obj *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::L1Obj >(nullptr);
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
      return GenerateInitInstanceLocal((::L1Obj*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1Obj*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_L1ObjColl(void *p = nullptr);
   static void *newArray_L1ObjColl(Long_t size, void *p);
   static void delete_L1ObjColl(void *p);
   static void deleteArray_L1ObjColl(void *p);
   static void destruct_L1ObjColl(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::L1ObjColl*)
   {
      ::L1ObjColl *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::L1ObjColl >(nullptr);
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
      return GenerateInitInstanceLocal((::L1ObjColl*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::L1ObjColl*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EventObj::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *EventObj::Class_Name()
{
   return "EventObj";
}

//______________________________________________________________________________
const char *EventObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int EventObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventObj*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GenObj::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *GenObj::Class_Name()
{
   return "GenObj";
}

//______________________________________________________________________________
const char *GenObj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int GenObj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenObj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenObj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObj*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr GenObjColl::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *GenObjColl::Class_Name()
{
   return "GenObjColl";
}

//______________________________________________________________________________
const char *GenObjColl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int GenObjColl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenObjColl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenObjColl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::GenObjColl*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr L1Obj::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *L1Obj::Class_Name()
{
   return "L1Obj";
}

//______________________________________________________________________________
const char *L1Obj::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int L1Obj::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1Obj::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1Obj::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1Obj*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr L1ObjColl::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *L1ObjColl::Class_Name()
{
   return "L1ObjColl";
}

//______________________________________________________________________________
const char *L1ObjColl::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int L1ObjColl::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1ObjColl::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1ObjColl::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::L1ObjColl*)nullptr)->GetClass(); }
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

namespace ROOT {
   static TClass *vectorlElongsPdoublegR_Dictionary();
   static void vectorlElongsPdoublegR_TClassManip(TClass*);
   static void *new_vectorlElongsPdoublegR(void *p = nullptr);
   static void *newArray_vectorlElongsPdoublegR(Long_t size, void *p);
   static void delete_vectorlElongsPdoublegR(void *p);
   static void deleteArray_vectorlElongsPdoublegR(void *p);
   static void destruct_vectorlElongsPdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<long double>*)
   {
      vector<long double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<long double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<long double>", -2, "vector", 389,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<long double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElongsPdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<long double>*)nullptr)->GetClass();
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
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
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
   static void *new_vectorlEboolgR(void *p = nullptr);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 596,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<bool>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bool>*)nullptr)->GetClass();
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
   static TClass *vectorlEL1ObjgR_Dictionary();
   static void vectorlEL1ObjgR_TClassManip(TClass*);
   static void *new_vectorlEL1ObjgR(void *p = nullptr);
   static void *newArray_vectorlEL1ObjgR(Long_t size, void *p);
   static void delete_vectorlEL1ObjgR(void *p);
   static void deleteArray_vectorlEL1ObjgR(void *p);
   static void destruct_vectorlEL1ObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<L1Obj>*)
   {
      vector<L1Obj> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<L1Obj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<L1Obj>", -2, "vector", 389,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<L1Obj>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEL1ObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<L1Obj>*)nullptr)->GetClass();
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
   static void *new_vectorlEGenObjgR(void *p = nullptr);
   static void *newArray_vectorlEGenObjgR(Long_t size, void *p);
   static void delete_vectorlEGenObjgR(void *p);
   static void deleteArray_vectorlEGenObjgR(void *p);
   static void destruct_vectorlEGenObjgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<GenObj>*)
   {
      vector<GenObj> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<GenObj>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<GenObj>", -2, "vector", 389,
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
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<GenObj>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEGenObjgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<GenObj>*)nullptr)->GetClass();
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
  void TriggerDictionaryInitialization_libOMTFObjects_Impl() {
    static const char* headers[] = {
"EventObj.h",
"GenObj.h",
"GenObjColl.h",
"L1Obj.h",
"L1ObjColl.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.26.04-edd28/x86_64-centos7-gcc11-opt/include",
"/scratch/Magisterka/Current_work/RootAnalysis/GMT/DataFormats/src/../include",
"/scratch/Magisterka/Current_work/RootAnalysis/GMT/DataFormats/src",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.26.04-edd28/x86_64-centos7-gcc11-opt/include/",
"/scratch/Magisterka/Current_work/RootAnalysis/build/GMT/DataFormats/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libOMTFObjects dictionary forward declarations' payload"
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
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libOMTFObjects dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "EventObj.h"
#include "GenObj.h"
#include "GenObjColl.h"
#include "L1Obj.h"
#include "L1ObjColl.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EventObj", payloadCode, "@",
"GenObj", payloadCode, "@",
"GenObjColl", payloadCode, "@",
"L1Obj", payloadCode, "@",
"L1ObjColl", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libOMTFObjects",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libOMTFObjects_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libOMTFObjects_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libOMTFObjects() {
  TriggerDictionaryInitialization_libOMTFObjects_Impl();
}
