#ifndef  __ColumnStencil_h
#define  __ColumnStencil_h

#include  "stencilinterface.h"
#include  "sparsestructureinterface.h"
#include  "gascoigne.h"


namespace Gascoigne
{

/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnStencil

////
////
/////////////////////////////////////////////

class ColumnStencil : public virtual StencilInterface
{
protected:

  IntVector   scol, sstart;

  void _RangeErrorStartStop(int i, const std::vector<int>& vec) const 
  {
    assert(i>=0);
    assert(i<=sstart.size());
  }  

public:

//
////  Con(De)structor 
//

  ColumnStencil() : StencilInterface() {}
  ~ColumnStencil() {}

  const IntVector&  col()    const { return scol; }
        IntVector&  col()          { return scol; }
  const IntVector&  start()  const { return sstart; }
        IntVector&  start()        { return sstart; }

  int  n()            const { return sstart.size()-1;}
  int  nentries()     const { return scol.size();}
  int  rowsize(int i) const { return sstart[i+1]-sstart[i];}

        int&  col(int pos)           { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
  const int&  col(int pos)     const { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
        int&  start(int i)           {  _RangeErrorStartStop(i,sstart); return sstart[i]; } 
  const int&  start(int i)     const {  _RangeErrorStartStop(i,sstart); return sstart[i]; } 
        int&  stop(int i)            {  _RangeErrorStartStop(i,sstart); return sstart[i+1]; } 
  const int&  stop(int i)      const {  _RangeErrorStartStop(i,sstart); return sstart[i+1]; } 

  void memory(const SparseStructureInterface*);
  void memory(int n, int nt);
  
  virtual int Find(int i, int j) const
    {
      for(int pos=start(i); pos<stop(i); pos++)
        {
          if (col(pos)==j) return pos;
        }
      std::cerr << "UnstructuredStencil::Find()";
      std::cerr << "no such coupling: "<< i <<" "<<j<<std::endl;
      abort();
      return -1;
    }

  std::ostream& Write(std::ostream& os) const;
};
}

#endif
