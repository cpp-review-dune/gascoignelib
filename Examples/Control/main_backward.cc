#include  "problemdescriptorbackward.h"
#include  "localtimeloop.h"

using namespace std;
using namespace Gascoigne;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("forward.param");

  if(argc!=5) {
    cerr << "usage: \"forward: name name first last\"";
  }
  string  iname(argv[1]);
  string  name(argv[2]);
  int first = atoi(argv[3]);
  int last = atoi(argv[4]);
  
  cerr << "************************************\n";
  cerr << "name first last:\t" << name << " " << first << " " << last << endl;
  cerr << "************************************\n";

  ProblemDescriptorBackward LPD;
  LPD.BasicInit(&paramfile);

  LocalTimeLoop loop;
  loop.BasicInit(&paramfile);
  loop.backward(iname,name,first,last,&LPD);

  return 0;
}
