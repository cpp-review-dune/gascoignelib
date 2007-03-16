/*----------------------------   faceq2.h     ---------------------------*/
/*      $Id$                 */
#ifndef __faceq2_H
#define __faceq2_H
/*----------------------------   faceq2.h     ---------------------------*/


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
    class FaceQ2 : public virtual FaceDiscretization
    {
    protected:

      void TransformationFace(FemInterface::Matrix& T1,FemInterface::Matrix& T2,int f) const;
	    
    public:
      FaceQ2() : FaceDiscretization() {}
      virtual ~FaceQ2() {}
      
      //
      //// Functions called from the Solver
      //
      virtual std::string GetName() const { return "Face Q2"; }
      
      virtual void BasicInit(const ParamFile* pf);
      virtual void ReInit   (const MeshInterface* M);
      virtual void build_faces();
    };
}



/*----------------------------   faceq2.h     ---------------------------*/
/* end of #ifndef __faceq2_H */
#endif
/*----------------------------   faceq2.h     ---------------------------*/
