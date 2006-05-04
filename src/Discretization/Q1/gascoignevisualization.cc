#include  "gascoignevisualization.h"
#include  "componentinformation.h"
#include  "compose_name.h"


/*-----------------------------------------*/

namespace Gascoigne
{

void GascoigneVisualization::AddVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  assert(mesh);
  assert(v);
  _v=v;

  int ncomp = v->ncomp();
  assert(ncomp);

  VDI.Clear();

  VD .SetGlobalVector(v);

  //needed for when generating vectors
  CI->SetDimension(mesh->dimension());
  
  //VDI.AddScalars(ncomp);
  {
    int         ncomp2 = CI->GetNScalars();
    if(ncomp==ncomp2){
      std::string s_name;
      for(int i=0;i<ncomp;i++){
        CI->GetScalarName(i,s_name);
        VDI.AddScalar(i,s_name,i);
      }
    }else{
      std::string s_name;
      for(int i=0;i<ncomp;i++){
        s_name="u";
        compose_name_without_dot(s_name,i); 
        VDI.AddScalar(i,s_name,i);
      }
    }
  }  

  //add vectors
  {
    int             nvectors    = CI->GetNVectors();
    fixarray<3,int> fa_vectorindices;
    std::string     s_name;
    for(int i=0;i<nvectors;i++){
      CI->GetVectorName(i,s_name);
      CI->GetVectorIndices(i,fa_vectorindices);
      VDI.AddVector(i,s_name,fa_vectorindices);
    }
  }
}

/*-----------------------------------------*/

void GascoigneVisualization::AddVector(const GlobalVector* v) 
{
  assert(mesh);
  assert(v);
  _v=v;

  int ncomp = v->ncomp();
  assert(ncomp);

  VDI.Clear();

  VD .SetGlobalVector(v);
  VDI.AddScalars(ncomp);

  int vector_index=0;
  if (ncomp>2) 
    {
      fixarray<3,int> ff;
      int dim = mesh->dimension();
      if (dim==3)
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = 3;
	}
      else if (dim==2)
	{
	  ff[0] = 1; ff[1] = 2; ff[2] = -1;
	}
      VDI.AddVector(vector_index,"v",ff);
      vector_index++;
    }
}

/*-----------------------------------------*/

void GascoigneVisualization::AddPointVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  AddVector(CI,v);

  SetPointData(&VD);
  SetPointDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddPointVector(const GlobalVector* v) 
{
  AddVector(v);

  SetPointData(&VD);
  SetPointDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddCellVector(const ComponentInformation* CI, const GlobalVector* v) 
{
  AddVector(CI,v);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);
}

/* -------------------------------------------------------*/

void GascoigneVisualization::AddCellVector(const GlobalVector* v) 
{
  AddVector(v);

  SetCellData(&VD);
  SetCellDataInfo(&VDI);
}
}
