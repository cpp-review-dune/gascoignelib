#ifndef __visudatainfo_h
#define __visudatainfo_h

#include  "fixarray.h"
#include  "visudata.h"

#include  <map>
#include  <string>

/*-------------------------------------------------------------------------*/

class VisuDataInfo
{
 protected:

  std::map<std::string,int>              scalars;
  std::map<std::string,fixarray<3,int> > vectors;

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
    scalars.clear();
    vectors.clear();
  }


  void AddScalar(const std::string& name, int i)                    {scalars[name]=i;}
  void AddVector(const std::string& name, const fixarray<3,int>& i) {vectors[name]=i;}

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

#endif
