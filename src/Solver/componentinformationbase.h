#ifndef  __ComponentInformationBase_h
#define  __ComponentInformationBase_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"

#include  "componentinformation.h"


#define CLASS ComponentInformationBase
namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ComponentInformationBase

  ///
  ///
  /////////////////////////////////////////////

  class CLASS : public ComponentInformation
  {
    private:
      
    protected:

    public:
      CLASS():ComponentInformation() {}
      virtual ~CLASS() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const;


      virtual const int GetNScalars     () const;
      virtual void      GetScalarName   (int i, std::string& s_name) const;
      virtual const int GetNVectors     () const;
      virtual void      GetVectorName   (int i, std::string& s_name) const;
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const;
  };
}

#undef CLASS // ComponentInformationBase

#endif // __ComponentInformationBase_h
