/*----------------------------   dgequation.h     ---------------------------*/
/*      $Id:$                 */
#ifndef __dgequation_H
#define __dgequation_H
/*----------------------------   dgequation.h     ---------------------------*/


#include "equation.h"

namespace Gascoigne
{

  class DGEquation : virtual public Equation
  {
  public:

    virtual void SetFemDataMaster(FemData& q) const {  }
    virtual void SetFemDataSlave(FemData& q) const {  }
    

    virtual void point_edge(bool internaledge,
			    double h,
                            const FemFunction &U1,
                            const FemFunction &U2,
                            const Vertex2d &v,
                            const Vertex2d &n) const
    {
    }
    virtual void point_edge(bool internaledge,
			    double h,
                            const FemFunction &U1,
                            const FemFunction &U2,
                            const Vertex3d &v,
                            const Vertex3d &n) const
    {
    }


    virtual void EdgeForm1(VectorIterator b,
			   const FemFunction& U1,
			   const FemFunction& U2,
			   const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeForm not written" << std::endl;
      abort();
    }
    virtual void EdgeForm2(VectorIterator b,
			   const FemFunction& U1,
			   const FemFunction& U2,
			   const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeForm not written" << std::endl;
      abort();
    }
    virtual void EdgeMatrix11(EntryMatrix& A,
			      const FemFunction& U1,
			      const FemFunction& U2,
			      const TestFunction& M,
			      const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeMatrix not written" << std::endl;
      abort();
    }
    virtual void EdgeMatrix12(EntryMatrix& A,
			      const FemFunction& U1,
			      const FemFunction& U2,
			      const TestFunction& M,
			      const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeMatrix not written" << std::endl;
      abort();
    }
    virtual void EdgeMatrix21(EntryMatrix& A,
			      const FemFunction& U1,
			      const FemFunction& U2,
			      const TestFunction& M,
			      const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeMatrix not written" << std::endl;
      abort();
    }
    virtual void EdgeMatrix22(EntryMatrix& A,
			      const FemFunction& U1,
			      const FemFunction& U2,
			      const TestFunction& M,
			      const TestFunction& N) const
    {
      std::cerr << "DGEquation::EdgeMatrix not written" << std::endl;
      abort();
    }
    
  };
}


/*----------------------------   dgequation.h     ---------------------------*/
/* end of #ifndef __dgequation_H */
#endif
/*----------------------------   dgequation.h     ---------------------------*/
