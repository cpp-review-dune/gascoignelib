#ifndef  __Q1Lps3d_h
#define  __Q1Lps3d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Lps3d

////
////
/////////////////////////////////////////////

#include  "q13d.h"
#include  "q1lpsstab.h"

namespace Gascoigne
{

/*----------------------------------------------*/

class Q1Lps3d : public Q13d
{
 protected:

  Q1LpsStab3d   S;

public:

//
////  Con(De)structor 
//

  Q1Lps3d() : Q13d() {}
  ~Q1Lps3d() {}

  std::string GetName() const {return "Q1Lps3d";}
  
  void BasicInit(const ParamFile* paramfile);
  void ReInit   (const MeshInterface* M);
  void Structure(SparseStructureInterface* SI) const;
  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AddNodeVector(const std::string& name, const GlobalVector* q) const 
    {
      Q13d::AddNodeVector(name,q);
      S.    AddNodeVector(name,q);
    }
  void DeleteNodeVector(const std::string& name) const 
    {
      Q13d::DeleteNodeVector(name);
      S.    DeleteNodeVector(name);
    }
  void AddCellVector(const std::string& name, const GlobalCellVector* q) const 
    {
      Q13d::AddCellVector(name,q);
      S.    AddCellVector(name,q);
    }
  void DeleteCellVector(const std::string& name) const
    { 
      Q13d::DeleteCellVector(name);
      S.    DeleteCellVector(name);
    }
  void AddParameterVector(const std::string& name, const GlobalParameterVector* q) const 
    {
      Q13d::AddParameterVector(name,q);
      S.    AddParameterVector(name,q);
   }
  void DeleteParameterVector(const std::string& name) const
    {
      Q13d::DeleteParameterVector(name);
      S.    DeleteParameterVector(name);
    }
};

}
#endif