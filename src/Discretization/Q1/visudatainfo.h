#ifndef __visudatainfo_h
#define __visudatainfo_h

#include  "fixarray.h"
#include  "visudata.h"

#include  <map>
#include  <string>

/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
class VisuDataInfo
{
 protected:

  std::map<std::string,int>              scalars;
  std::map<std::string,fixarray<3,int> > vectors;
  std::map<std::string,int>              scalar_order;
  std::map<std::string,int>              vector_order;

 public:

  typedef std::map<std::string,int>::const_iterator                siterator;
  typedef std::map<std::string,fixarray<3,int> >::const_iterator   viterator;

  VisuDataInfo() {}
  VisuDataInfo(int ncomp) { AddScalars(ncomp);}
  VisuDataInfo(const VisuData& V, std::string def="U");
  VisuDataInfo(const VisuDataInfo& V) : scalars(V.Scalars()), vectors(V.Vectors()) {}
  VisuDataInfo& operator=(const VisuDataInfo& V);

  bool operator!=(const VisuDataInfo& V) const;

  void Clear() {
    scalar_order.clear();
    vector_order.clear();
    scalars.clear();
    vectors.clear();
  }

  siterator GetSIterator(int i) { 
    for(siterator p = sbegin() ; p!= send() ; p++){
      std::string s = p->first;
      if ( scalar_order[s]==i ) return p;
    }
    assert(0);
  }
  viterator GetVIterator(int i) {
    for(viterator p = vbegin() ; p!= vend() ; p++){
      std::string s = p->first;
      if ( vector_order[s]==i ) return p;
    }
    assert(0);
  }

  void AddScalar(int index,const std::string& name, int i)                    {scalar_order[name]=index;scalars[name]=i;}
  void AddVector(int index,const std::string& name, const fixarray<3,int>& i) {vector_order[name]=index;vectors[name]=i;}

  void AddScalars(int ncomp, std::string def="U");

  int nscalars() const {return scalars.size();}
  int nvectors() const {return vectors.size();}

  const std::map<std::string,int>&              Scalars() const {return scalars;}
  const std::map<std::string,fixarray<3,int> >& Vectors() const {return vectors;}

  siterator sbegin() const {return scalars.begin();}
  siterator send  () const {return scalars.end();}
  viterator vbegin() const {return vectors.begin();}
  viterator vend  () const {return vectors.end();}
};
}

#endif
