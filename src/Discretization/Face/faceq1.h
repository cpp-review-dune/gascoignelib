/*----------------------------   faceq1.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceq1_H
#define __faceq1_H
/*----------------------------   faceq1.h     ---------------------------*/



#include  "facediscretization.h"

namespace Gascoigne
{
  
  /////////////////////////////////////////////
  ///
  ///@brief
  ///  ... comments DiscretizationInterface

  ///
  ///
  /////////////////////////////////////////////

  template<int DIM>
    class FaceQ1 : public virtual FaceDiscretization
    {
    protected:

      void TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const;
      
    public:
      FaceQ1() : FaceDiscretization() {}
      virtual ~FaceQ1() {}
      
      //
      //// Functions called from the Solver
      //
      virtual std::string GetName() const { return "Face Q1"; }
      
      virtual void BasicInit(const ParamFile* pf);
      virtual void ReInit   (const MeshInterface* M);
      virtual void build_faces();
    };
}




/*----------------------------   faceq1.h     ---------------------------*/
/* end of #ifndef __faceq1_H */
#endif
/*----------------------------   faceq1.h     ---------------------------*/
