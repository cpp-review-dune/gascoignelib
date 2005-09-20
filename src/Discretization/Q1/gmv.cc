#include  "visualization.h"
#include  "errormacros.h"

using namespace std;


/********************************************************************/

namespace Gascoigne
{
void Visualization::gmv(const string& bname) const
{
  string name = bname;
  name += ".gmv";

  ofstream file(name.c_str());
  FILE_ERROR(file,name);
  
  file << "gmvinput ascii" << endl;
  file << "nodes " << mesh->nnodes() << endl;
  
  output_vertexs_by_component(file);

  file << endl << "cells " << mesh->ncells() << endl;

  output_quads(file,"quad 4 ");
  output_hexs (file,"hex 8 ");

  if (PointData)
    {
      CheckPointData();
      file << "variable" << endl;
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
	{
	  file << p->first <<" 1" << endl;
	  output_solution(file,p->second);
	  file << endl;
	}
     file << "endvars" << endl;
    }
  file << "endgmv" << endl;
  file.close();
 if (showoutput) cout << "[" << name << "]\n";

 if(compress)
 {
   string command = "gzip " + name;
   system(command.c_str());
 }
}
}
