#include  "problemdescriptorforward.h"
#include  "localtimeloop.h"

using namespace std;

/*---------------------------------------------------*/

int main(int argc, char** argv)
{
  ParamFile paramfile("forward.param");

  if(argc!=4) {
    cerr << "usage: \"forward: name first last\"";
  }
  string  name(argv[1]);
  int first = atoi(argv[2]);
  int last = atoi(argv[3]);
  
  cerr << "************************************\n";
  cerr << "name first last:\t" << name << " " << first << " " << last << endl;
  cerr << "************************************\n";

  ProblemDescriptorForward LPD;
  LPD.BasicInit(&paramfile);

  LocalTimeLoop loop;
  loop.BasicInit(&paramfile);
  loop.forward(name,first,last,&LPD);

  return 0;
}
