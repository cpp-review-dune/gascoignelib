#ifndef  __ColumnStencil_h
#define  __ColumnStencil_h


/////////////////////////////////////////////
////
////@brief
////  ... comments ColumnStencil

////
////
/////////////////////////////////////////////

#include  "stencilinterface.h"
#include  "sparsestructureinterface.h"
#include  "nvector.h"

class ColumnStencil : public virtual StencilInterface
{
private:


protected:

  void _RangeErrorStartStop(int i, const std::vector<int>& vec) const {
/*     if( !((i>=0)&&(i+1<sstart.size())) ) */
/*       { */
/* 	cerr << "start/stop out of range: i="<<i<< " vec.size() "<<vec.size()<<endl; */
/* 	assert(0); */
/*       }  */
  }  
  nvector<int>   scol, sstart;

public:


//
////  Con(De)structor 
//

  ColumnStencil() : StencilInterface() {}
  ~ColumnStencil() {}

  const nvector<int>&  col()    const { return scol; }
        nvector<int>&  col()          { return scol; }
  const nvector<int>&  start()  const { return sstart; }
        nvector<int>&  start()        { return sstart; }

  int  n()        const { return sstart.size()-1;}
  int  nentries() const { return scol.size();}
  int  rowsize(int i) const { return sstart[i+1]-sstart[i];}

        int&  col(int pos)           { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
  const int&  col(int pos)     const { assert((pos>=0)&&(pos<scol.size())); return scol[pos]; } 
        int&  start(int i)           {  _RangeErrorStartStop(i,sstart); return sstart[i]; } 
  const int&  start(int i)     const { _RangeErrorStartStop(i,sstart); return sstart[i]; } 
  int&  stop(int i)                  {  _RangeErrorStartStop(i,sstart); return sstart[i+1]; } 
  const int&  stop(int i)      const {  _RangeErrorStartStop(i,sstart); return sstart[i+1]; } 

  void memory(int n, int nt);
  void memory(const SparseStructureInterface*);
  
  virtual int Find(int i, int j) const
    {
      bool error = 1;
      int pos;
      for(pos=start(i);pos<stop(i);pos++)
	{
	  if(col(pos)==j)   
	    {
	      error = 0;
	      break;
	    }
	}
      if(error) 
	{
	  std::cerr << "UnstructuredStencil::Find()";
	  std::cerr << "no such coupling: "<<i <<" "<<j<<std::endl;
	  abort();
	  return -1;
	}
      return pos;
    }

  std::ostream& Write(std::ostream& os) const
  {
    os << n() << "\t" << nentries()<<std::endl<<std::endl;
    os << sstart<<std::endl<<std::endl;
    for(int i=0;i<n();i++)
      {
	for(int pos=start(i);pos<stop(i);pos++)
	  {
	    os << col(pos) << " ";
	  }
	os << std::endl;
      }
    return os;
  }
};


#endif
