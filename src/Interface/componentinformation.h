#ifndef  __ComponentInformation_h
#define  __ComponentInformation_h

#include  "gostream.h"
#include  "stringutil.h"
#include  "filescanner.h"
#include  "application.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments ComponentInformation

  ///
  ///
  /////////////////////////////////////////////

  class ProblemDescriptorInterface;
  class ComponentInformation : public virtual Application
  {
    private:
      
    protected:
      mutable int _i_dimension;
      ProblemDescriptorInterface* PDI;

    public:
      ComponentInformation() { PDI=NULL;}
      virtual ~ComponentInformation() {}
  
      virtual void BasicInit(const ParamFile* pf) {}

      virtual std::string GetName() const=0;
      virtual int         GetDimension()                const { return _i_dimension;      };
      virtual void        SetDimension(int i_dimension) const { _i_dimension=i_dimension; };
      ProblemDescriptorInterface*& GetProblemDescriptorInterface()       { return PDI;};
      ProblemDescriptorInterface*  GetProblemDescriptorInterface() const { return PDI;};

      const int ncomp   () const { return GetNScalars(); };
      const int GetNcomp() const { return GetNScalars(); };

      virtual const int GetNScalars     () const=0;
      virtual void      GetScalarName   (int i, std::string& s_name) const=0;
      virtual const int GetNVectors     () const=0;
      virtual void      GetVectorName   (int i, std::string& s_name) const=0;
      virtual void      GetVectorIndices(int i, fixarray<3,int>& fa_vectorindices) const=0;
  };
}

#endif // __ComponentInformation_h
