#include  "visualization.h"
#include  "errormacros.h"

using namespace std;


/********************************************************************/

namespace Gascoigne
{
void Visualization::vu(const string& bname) const
{
  string name = bname;
  name += ".pie";

  string element;
  if (mesh->dimension()==2)
    {
      element = "LagrQuadr04";
    }
  else
    {
      element = "LagrHexae08";
    }
  
  ofstream file(name.c_str());
  FILE_ERROR(file,name);

  file << "MAILLAGE MonMaillage( ) =\n{" << endl;
  file << "   ZONE Zone( " << element << ", Nodes, Elements );\n};\n\n";

//   if (symmetry)
//     {
//       file << "IMAGE Mirror( Z=0, Graphe1Plein,Princ ) =\n";
//       file << "{\n  NCopies   2;\n  DRotation 180,1,0,0;\n};\n\n";
//     }
  file << "CHAMP Nodes( ) =\n{\n";
  
  output_vertexs(file);
  file << "};\n\n";

  file << "CHAMP<int> Elements( ) =\n{\n";
  if (mesh->dimension()==3)
    {
      output_hexs(file);
    }
  else
    {
      output_quads(file);
    }

  file << "};\n\n";

  if (PointData)
    {
      CheckPointData();
      file << "SOLUTION Solution( ) =\n{\n";
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
	{
	  file << " VARIABLE " << p->first << "( " << element;
	  file << ", " << p->first;
	  file << ", Elements, Zone );\n";
	}
      file << "};\n\n";
      for(VisuDataInfo::siterator p=PointDataInfo->sbegin();p!=PointDataInfo->send();++p)
	{
	  file << "CHAMP " << p->first << "( ) = {\n";
	  output_solution(file,p->second);
	  file << "};\n\n";
	}
    }
  file.close();

 if(compress)
 {
   string command = "gzip -f " + name;
   system(command.c_str());
 }
}
}
