#include  "visualization.h"
#include  "errormacros.h"


using namespace std;

/********************************************************************/

void Visualization::avs(const string& bname) const
{
  string name = bname;
  name += ".inp";
  
  ofstream out(name.c_str());
  FILE_ERROR(out,name);
  
  int nc = CheckPointData();
  
  bool CreateZero=0;
  for(VisuDataInfo::viterator p=PointDataInfo->vbegin();p!=PointDataInfo->vend();++p)
    {
      if(p->second[2]==-1)
	{
	  CreateZero=1;
	  continue;
	}
    }
  out << mesh->nnodes() << " " << mesh->ncells();
  out << " " << nc << " 0 0 " << endl;
  
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  out<< i+1 << " " << mesh->vertex2d(i) << " " << 0 << endl;
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  out<< i+1 << " " << mesh->vertex3d(i) << endl;
	}
    }
  if (mesh->dimension()==2)
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  out << c+1 << " 1 quad ";
	  for(int i=0;i<4;i++)
	    {
	      out << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  out << endl; 
	} 
    }
  else
    {
      for (int c=0; c<mesh->ncells(); c++)
	{
	  out << c+1 << " 1 hex ";
	  for(int i=0;i<8;i++)
	    {
	      out << mesh->vertex_of_cell(c,i)+1 << " "; 
	    }
	  out << endl; 
	}     
    }
  if (PointData)
    {
      out << nc;
      for (int c=0; c<nc; c++) {out << " 1";}
      out << endl;
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
       {
	 out << p->first << "," << endl;
       }
      if(CreateZero)  out << "null,";
	 
      for (int ind=0; ind<PointData->visun(); ind++)
	{
	  out << endl << ind;
	  for(VisuDataInfo::siterator p=PointDataInfo->sbegin();
	      p!=PointDataInfo->send();++p)
	    {
	      out << " " << PointData->visudata(ind,p->second);
	    }
	  if(CreateZero) out << " 0";
	}
      out << endl;
    }
  out.close();
 if (showoutput) cout << "[" << name << "]\n";
}
