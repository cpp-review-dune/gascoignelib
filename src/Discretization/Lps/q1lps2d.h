#ifndef  __Q1Lps2d_h
#define  __Q1Lps2d_h

/////////////////////////////////////////////
////
////@brief
////  ... comments Q1Lps2d

////
////
/////////////////////////////////////////////

#include  "q12d.h"
#include  "q1lpsstab.h"

namespace Gascoigne
{

/*----------------------------------------------*/

class Q1Lps2d : public Q12d
{
 protected:

  typedef Gascoigne::GlobalVector  GlobalVector; 

  Q1LpsStab2d   S;

public:

//
////  Con(De)structor 
//

  Q1Lps2d() : Q12d() {}
  ~Q1Lps2d() {}

  std::string GetName() const {return "Q1Lps2d";}
  
  void BasicInit(const Gascoigne::ParamFile* paramfile);
  void ReInit   (const MeshInterface* M);
  void Structure(SparseStructureInterface* SI) const;
  void Form(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void Matrix(MatrixInterface& A, const GlobalVector& u, const Equation& EQ, double d) const;
  void StabForm(GlobalVector& f, const GlobalVector& u, const Equation& EQ, double d) const;
  void AddNodeVector(const std::string& name, const Gascoigne::GlobalVector* q) const 
    {
      Q12d::AddNodeVector(name,q);
      S.    AddNodeVector(name,q);
    }
  void DeleteNodeVector(const std::string& name) const 
    {
      Q12d::DeleteNodeVector(name);
      S.    DeleteNodeVector(name);
    }
};

}

#endif
