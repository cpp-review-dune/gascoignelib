#ifndef  __LocalTimeLoop_h
#define  __LocalTimeLoop_h

#include  "stdtimeloop.h"


/////////////////////////////////////////////
////
////@brief
////  ... comments LocalTimeLoop

////
////
/////////////////////////////////////////////



class LocalTimeLoop : public StdTimeLoop
{
private:

  vector<GlobalVector>  dat_node;

protected:


public:


//
////  Con(De)structor 
//

  LocalTimeLoop() : StdTimeLoop() {}
  void BasicInit(const ParamFile* paramfile);

  void AddNodeVector(string filename);
  void NewMesh(const ProblemDescriptorInterface* PD);

  void init(string,int);
  void forward(string iname, int first, int last, const ProblemDescriptorInterface* PD);
  void backward(string iname, string name, int first, int last, const ProblemDescriptorInterface* PD);

};


#endif
