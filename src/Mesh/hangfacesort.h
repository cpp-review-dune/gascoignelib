#ifndef __hangfacesort_h
#define __hangfacesort_h

#include "facemanager.h"

/*---------------------------------------------------*/

namespace Gascoigne
{
class HangFaceSort{

protected:

  const FaceManager& HR;

public:

  HangFaceSort(const FaceManager& H) : HR(H) {}
  bool operator() (int i, int j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};

/*---------------------------------------------------*/

class HangFaceSort2{

protected:

  const FaceManager& HR;

public:

  HangFaceSort2(const FaceManager& H) : HR(H) {}
  bool operator() (const Edge& i, const Edge& j) const
    {
      return !HR.EdgeIsHanging(i) && HR.EdgeIsHanging(j);
    }
};
}

#endif
