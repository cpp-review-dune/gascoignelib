#include  "starter.h"
#include  "assert.h"

using namespace std;

//-----------------------------------------------------------------//

Starter::Starter(int argc, char** argv, const std::string& paramfile) 
{    
  _paramfile = new ParamFile;
  if(argc>=2) 
    {
      _erase = 0;
      _paramfile->SetName(argv[1]);
    } 
  else 
    {
      string pname;
      _erase = 1;
      pname = "_" + paramfile;
      _paramfile->SetName(pname);
      //std::cout << "Paramfile ist: " << pname << std::endl;
      {
	string cmd("rm -f ");
	cmd += pname;
	if( system(cmd.c_str())) exit(-1);
      }
      {
	string cmd("cp ");
	cmd += paramfile;
	cmd += " ";
	cmd += pname;
	if( system(cmd.c_str())) exit(-1);
      }
      {
	if( !system("ls common.param") ) 
	  {
	    string cmd("cat ");
	    cmd += "common.param";
	    cmd += " >> ";
	    cmd += pname;
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
      assert(_paramfile);
      cmd += _paramfile->GetName();
      if( system(cmd.c_str())) exit(-1);
    }
}


