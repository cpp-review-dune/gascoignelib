#ifndef __visu_eps_h
#define __visu_eps_h

#include "patchmesh.h"
#include "string"
#include <map>
#include  "gascoigne.h"

/*-------------------------------------------------------------------------*/

namespace Gascoigne
{
class VisuEPS
{
  protected:  

  typedef std::pair<int,int> Line;
  typedef nvector<IntSet> Lines;

  const PatchMesh* M;

  Lines lines;
  int   n_lines;
  Vertex2d offset;
  
  // EPS optionen
  std::map<int,int>    INTOPT;
  std::map<int,double> DOUBLEOPT;

  // sort line p (left,bottom) first
  void Lexiko(Line& p) const;

  // test if vertices a,b,c are aligned straightly
  bool InLine(int a,int b,int c) const;
  
  void CombineLines();
    
  public:

  VisuEPS();

  /**
   * Options for output:
   *
   * WRITE_PATCH:
   *   0 : write cells  (default)
   *   1 : write patchs
   *
   * LINEWIDTH:
   *   width of lines in pt. (0.1)
   *
   * WIDTH:
   *   horizontal size of output (int pt) (300)
   *
   * COMBINE_LINES
   *   1: straightly aligned lines are combined to one (def)
   *
   **/
  enum EPSOptions { WRITE_PATCH, LINEWIDTH, WIDTH, COMBINE_LINES };
  
  void SetOption(EPSOptions o, int v);
  void SetOption(EPSOptions o, double v);
  
  void SetMesh(const MeshInterface& PM)  { 
    const PatchMesh* PMP = dynamic_cast<const PatchMesh*>(&PM);
    assert(PMP);
    M = PMP;}
  void WriteGrid(std::string fname, int iter);
};
}

#endif
