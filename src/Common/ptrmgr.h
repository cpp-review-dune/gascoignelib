#ifndef  __PtrMgr_h
#define  __PtrMgr_h


/////////////////////////////////////////////
////
////@brief
////  ... comments PtrMgr

////
////
/////////////////////////////////////////////


#include  <string>
#include  <map>

using namespace std;

template<class T>
class PtrMgr
{
public:

  typedef typename map<string,T*>::const_iterator const_iterator;
  typedef typename map<string,T*>::iterator iterator;

private:

  map<string,T*>  _m;

protected:


public:


//
////  Con(De)structor 
//
PtrMgr() {}
~PtrMgr() {
  for(iterator p=_m.begin();p!=_m.end();p++) {
    if(p->second==NULL) {delete p->second; p->second=NULL;}
  }
}

 const T* Get(const string& name) const {
   const_iterator p = _m.find(name);
   assert(p!=_m.end());
   return *p;
 }
 T*& GetPointer(const string& name) {
   return _m[name];
 }
 const_iterator first() const {
   return _m.begin();
 }
 const_iterator last() const {
   return _m.end();
 }
};


#endif
