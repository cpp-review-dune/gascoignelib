#include  "visualization.h"
#include  "errormacros.h"


/********************************************************************/

void Visualization::avs(const std::string& bname) const
{
  std::string name = bname;
  name += ".inp";
  
  std::ofstream out(name.c_str());
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
  out << " " << nc << " 0 0 " << std::endl;
  
  if (mesh->dimension()==2)
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  out<< i+1 << " " << mesh->vertex2d(i) << " " << 0 << std::endl;
	}
    }
  else
    {
      for (int i=0; i<mesh->nnodes(); i++)
	{
	  out<< i+1 << " " << mesh->vertex3d(i) << std::endl;
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
	  out << std::endl; 
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
	  out << std::endl; 
	}     
    }
  if (PointData)
    {
      out << nc;
      for (int c=0; c<nc; c++) {out << " 1";}
      out << std::endl;
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
       {
	 out << p->first << "," << std::endl;
       }
      if(CreateZero)  out << "null,";
	 
      for (int ind=0; ind<PointData->visun(); ind++)
	{
	  out << std::endl << ind;
	  for(VisuDataInfo::siterator p=PointDataInfo->sbegin();
	      p!=PointDataInfo->send();++p)
	    {
	      out << " " << PointData->visudata(ind,p->second);
	    }
	  if(CreateZero) out << " 0";
	}
      out << std::endl;
    }
  out.close();
 if (showoutput) std::cout << "[" << name << "]\n";
}
