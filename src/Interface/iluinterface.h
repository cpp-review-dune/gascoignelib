#ifndef __iluinterface_h
#define __iluinterface_h

#include  "matrixinterface.h"

/*-------------------------------------------------------------*/

class IluInterface
{
  public:

  virtual ~IluInterface() {};

  virtual std::string GetName() const  { assert(0);}

  virtual void modify(int c, double s) { assert(0);}  
  virtual void zero()                  { assert(0);}                      
  virtual void compute_ilu ()          { assert(0);}           
  virtual void ReInit   (const SparseStructureInterface* A) { assert(0);}       
  virtual void ConstructStructure(const nvector<int>& perm, const MatrixInterface& A)
    { assert(0); }
  virtual void copy_entries(const MatrixInterface* A)      { assert(0); }  
  virtual void solve       (CompVector<double>& x) const   { assert(0); }        
  virtual void solve_transpose(CompVector<double>& x) const   { assert(0); }    
  virtual std::ostream& Write(std::ostream &s) const       { assert(0); }
};

#endif
