#include  "visudatainfo.h"
#include  "compose_name.h"

using namespace std;

/*-------------------------------------------------------------------------*/

VisuDataInfo& VisuDataInfo::operator=(const VisuDataInfo& V)
{
  scalars = V.Scalars();
  vectors = V.Vectors();
  return *this;
}

/*-------------------------------------------------------------------------*/

bool VisuDataInfo::operator!=(const VisuDataInfo& V) const
{
  return (scalars!=V.Scalars())||(vectors!=V.Vectors());
}

/*-------------------------------------------------------------------------*/

VisuDataInfo::VisuDataInfo(const VisuData& D, string def)
{
  for(int c=0;c<D.visucomp();c++)
    {
      string name(def);
      compose_name_without_dot(name,c);
      AddScalar(name,c);
    }
}

/*-------------------------------------------------------------------------*/

VisuDataInfo::VisuDataInfo(int nc, string def)
{
  for(int c=0;c<nc;c++)
    {
      string name(def);
      compose_name_without_dot(name,c);
      AddScalar(name,c);
    }
}
