#include  "monitor.h"
#include  <string>
#include  <fstream>
#include  <stdio.h>
#include  "filescanner.h"
#include  "stlio.h"

/*****************************************************************/

Monitor::Monitor() : prec(6)
{
  aos = "Gascoigne";
  bos = "Adaptive";
  directory = ".";
  ps  = 1;
}

/*******************************************************************/

void Monitor::PrintInfoSummary(const CGInfo& info) const
{
  std::cout.precision(3);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << info.statistics().totaliter() << "\t";

  std::cout.setf(std::ios::scientific,std::ios::floatfield);
  std::cout.precision(2);
  std::cout << " { ";
  std::cout << info.control().residual() << " ";
  std::cout << info.control().firstresidual() << " (";
  std::cout << info.control().correction() << ") ";
  std::cout << "}\t";

  std::cout.precision(2);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << info.statistics().rate() << "\n";
}

/*******************************************************************/

void Monitor::PrintInfoSummary(const NLInfo& nlinfo) const
{
  const CGInfo& cginfo = nlinfo.GetLinearInfo();

  std::cout.precision(3);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << nlinfo.statistics().totaliter() << " [";
  std::cout << cginfo.statistics().totaliter() << "] M ";
  std::cout << nlinfo.statistics().totalmatrix() << " ";

  std::cout.setf(std::ios::scientific,std::ios::floatfield);
  std::cout.precision(2);
  std::cout << " { ";
  std::cout << nlinfo.control().residual() << " ";
  std::cout << nlinfo.control().firstresidual() << " (";
  std::cout << nlinfo.control().correction() << ") ";
  std::cout << "}\t";

  std::cout.precision(2);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  std::cout << nlinfo.statistics().rate() << " ";
  std::cout << cginfo.statistics().rate() << "\n";
}

/*****************************************************************/

void Monitor::set_directory(const std::string& dir)
{
  //std::cout << "Monitor::set_directory()\tdirectory = " << dir << std::endl;
  directory = dir;
  texfile   = directory;
  texfile += "/TexResults";
  numfile   = directory;
  numfile += "/NumResults";
  protokoll = directory;
  protokoll += "/Protokoll";

  std::string command("mkdir -p ");
  command += dir;
  system(command.c_str());

  std::ofstream tfile( texfile.c_str() );
  if(!tfile)
    {
      error_io(texfile);
    }
  else
    {
      PrintHeader(tfile);
    }
  tfile.close();
  
  std::ofstream pfile( protokoll.c_str() );
  if(!pfile)
    {
      error_io(protokoll);
    }
  pfile.close();
  
  std::ofstream nfile( numfile.c_str() );
  if(!nfile)
    {
      error_io(numfile);
    }
  nfile.close();
  sprintf (message,"");

  //cout << "set_directory:" << protokoll <<  "#" << endl;
}

/*****************************************************************/

void Monitor::init(const ParamFile* pf, int c)
{
  control   = c;

  std::string dir;
  DataFormatHandler   DFH;
  DFH.insert("directory"  , &dir,  ".");
  DFH.insert("format"  , &format,  "latex");
  DFH.insert("header"  , &header);
  FileScanner FS(DFH, pf, "Monitor");

  if (header.size()>0) std::cout << "header:\n" << header << std::endl;

  if(dir!=".")  set_directory(dir);
}

/*****************************************************************/

void Monitor::error_io(const std::string& s) const
{
  std::cout << "Monitor::error_io()\n";
  std::cout << "Fehler beim Oeffnen von " << s << std::endl;
  std::cout << "directory = " << directory << std::endl;
  abort();
}

/*****************************************************************/

void Monitor::pre_monitor(char* s)
{
  if(control>2)
    {
      m_length[m_counter++] = strlen(s);
      std::cout << s;
    }
}

/*****************************************************************/

void Monitor::post_monitor()
{
  if(control>2)
    {
      int l = m_length[--m_counter];
      for(int i=0; i<l; i++) std::cout << "\b";
      for(int i=0; i<l; i++) std::cout << " " ;
      for(int i=0; i<l; i++) std::cout << "\b";
    }
}

/*****************************************************************/

void Monitor::print_message()
{
  std::ofstream pfile( protokoll.c_str(),std::ios::app);

  if(!pfile)
  {
    error_io(protokoll);
  }
  if (control) std::cout << message << std::endl;

  pfile << message << std::endl;

  sprintf (message,"");
  pfile.close();
}

/*******************************************************************/

void Monitor::failed_step()
{
  print_message();
  sprintf(message,"%s repeat",message);
  print_message();
}
  
/*******************************************************************/

void Monitor::mesh(int nc, int nl)
{
  sprintf(message,"%s (N=%5d[%2d])", message, nc, nl);
  print_message();
}

/*******************************************************************/

void Monitor::post_nonlinear(const dvector& Ju, double eta, int ncells,
			     int cg, int dcg)
{
  sprintf (message,"cg dnl: %d %d",cg,dcg);
  sprintf (message,"%s [J(u)=%10.7e eta=%7.3e]",message,Ju[0],eta);
  
  print_message();
}

/*******************************************************************/

void Monitor::pre_nonlinear(int i)
{
  print_message();
  sprintf(message,"----------------------------------------");
  print_message();
  sprintf(message,"%s-%s: %4d",aos.c_str(),bos.c_str(),i);
}

/*******************************************************************/

/*******************************************************************/

void Monitor::matrix_info(int i)
{
  if(!i)
    {
      for(int j=0; j<ps; j++) sprintf(message,"%s ",message);
    }
  else
    {
      if (_newmatrix) sprintf(message,"%sM",message);
      else            sprintf(message,"%s ",message);
    }
}

/*******************************************************************/

void Monitor::nonlinear_step(const CGInfo& cginfo, const NLInfo& nlinfo)
{
  int i = nlinfo.control().iteration();

  matrix_info(i);

  if (!ps || (i%ps)) return;

  int  iter = cginfo.control().iteration(); 

  double r = nlinfo.control().residual();
  sprintf(message,"%s%3d: %8.2e",message,i,r);
  if (i && (control>1))
    {
      double c = nlinfo.control().correction();
      sprintf(message,"%s %8.2e",message,c);
    }
  if (i)
    {
      double rr = nlinfo.statistics().rate();
      double lr = nlinfo.statistics().lastrate();
      int    p  = nlinfo.control().relax();
      
      sprintf(message,"%s [%4.2f %4.2f]",message,lr,rr);
      if (p)
	{
	  sprintf(message,"%s %1d",message,p);
	}
      else
	{
	  sprintf(message,"%s  ",message);
	}
      
      rr = cginfo.statistics().rate();
      lr = cginfo.control().residual();
      double lc = cginfo.control().correction();
      
      sprintf(message,"%s - %8.2e",message,lr);
      if (control>1)
	{
	  sprintf(message,"%s %8.2e",message,lc);
	}
      sprintf(message,"%s [%5.3f] {%2d",message,rr,iter);
      if (cginfo.control().status()=="running") sprintf(message,"%s*",message);
      sprintf(message,"%s}",message);
    }
  if (cginfo.control().status()=="diverged")
    {
      sprintf(message,"%s @",message);
    }
  print_message();
}

/*******************************************************************/

void Monitor::PrintResults(const std::string& s) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);
  if(format=="latex")   PrintAscii(file,s);
  // ausser betrieb
//   if(format=="html")    PrintHtml (file,iv,dv,"\n");
//   if(format=="ascii")   PrintHtml (file,iv,dv,"\n");
  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(int i) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);
  if(format=="latex")   PrintAscii(file,i);
  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(double d) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);
  if(format=="latex")   PrintAscii(file,d);
  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(const ivector& iv) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);
  if(format=="latex")   PrintAscii(file,iv);
  file.close();
}

/*******************************************************************/

void Monitor::PrintResults(const dvector& dv) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);
  if(format=="latex")   PrintAscii(file,dv);
  file.close();
}

/*******************************************************************/

void Monitor::Print(const std::string& s, std::string se) const
{
  std::ofstream   file(texfile.c_str(),std::ios::app);
  if(!file)  error_io(texfile);

  file << s << se;
  std::cout << s << se;
  file.close();
}

/*******************************************************************/

void Monitor::PrintHeader(std::ostream& os) const
{
  if(header.size()==0) return;
  for (int i=0; i<header.size()-1; i++) os << header[i] << " , "; 
  os << header[header.size()-1] << std::endl;
}


/*******************************************************************/

void Monitor::PrintAscii
(std::ostream& os, const ivector& iv) const
{
  os.precision(prec);
  os << iv;
}

/*******************************************************************/

void Monitor::PrintAscii
(std::ostream& os, const dvector& dv) const
{
  os.precision(prec);
  os.setf(std::ios::scientific,std::ios::floatfield);
  os << dv;
}

/*******************************************************************/

void Monitor::PrintAscii
(std::ostream& os, const std::string& s) const
{
  os << s;
}
