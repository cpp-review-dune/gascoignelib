#include  "visualization.h"
#include  "errormacros.h"


/********************************************************************/

void Visualization::gmv(const std::string& bname) const
{
  std::string name = bname;
  name += ".gmv";

  std::ofstream file(name.c_str());
  FILE_ERROR(file,name);
  
  file << "gmvinput ascii" << std::endl;
  file << "nodes " << mesh->nnodes() << std::endl;
  
  output_vertexs_by_component(file);

  file << std::endl << "cells " << mesh->ncells() << std::endl;

  output_quads(file,"quad 4 ");
  output_hexs (file,"hex 8 ");

  if (PointData)
    {
      CheckPointData();
      file << "variable" << std::endl;
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
	{
	  file << p->first <<" 1" << std::endl;
	  output_solution(file,p->second);
	  file << std::endl;
	}
     file << "endvars" << std::endl;
    }
  file << "endgmv" << std::endl;
  file.close();
 if (showoutput) std::cout << "[" << name << "]\n";
}
