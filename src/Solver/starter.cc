#include "starter.h"

using namespace std;

Starter::Starter(int argc, char** argv, const std::string& paramfile) 
{    
  if(argc>=2) 
    {
      _erase = 0;
      _paramfile = argv[1];
    } 
  else 
    {
      _erase = 1;
      _paramfile = "_" + paramfile;
      //std::cout << "Paramfile ist: " << _paramfile << std::endl;
      {
	string cmd("rm -f ");
	cmd += _paramfile;
	if( system(cmd.c_str())) exit(-1);
      }
      {
	string cmd("cp ");
	cmd += paramfile;
	cmd += " ";
	cmd += _paramfile;
	if( system(cmd.c_str())) exit(-1);
      }
      {
	if( !system("ls common.param") ) 
	  {
	    string cmd("cat ");
	    cmd += "common.param";
	    cmd += " >> ";
	    cmd += _paramfile;
	    if( system(cmd.c_str())) exit(-1);
	  }
	else
	  {
	    //	      cerr << "++++++++ no file \"common.param\" ++++++++\n";
	  }
      }
    }
}

/*-----------------------------------------*/

Starter::~Starter() 
{
  if (_erase==1)
    {
      string cmd("rm -f ");
      cmd += _paramfile;
      if( system(cmd.c_str())) exit(-1);
    }
}


