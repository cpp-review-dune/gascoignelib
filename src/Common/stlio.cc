#include  "stlio.h"
#include  "vertex.h"
#include  "compvector.h"


// /*-------------------------------------------------------------------*/

// std::ostream& operator<<(std::ostream &s, const std::map<int,fixarray<4,int> >& A)
// {
//   for(std::map<int,fixarray<4,int> >::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t"<<p->first<<"\t --> "<<p->second<<std::endl;
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// std::ostream& operator<<(std::ostream &s, const std::map<std::pair<std::string,std::string>,int>& A)
// {
//   for(std::map<std::pair<std::string,std::string>,int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first.first<<"\t"<<p->first.second<<"\t --> "<<p->second<<std::endl;
//     }
//   return s;
// }

// std::ostream& operator<<(std::ostream &s, const std::map<std::string,int>& A)
// {
//   for(std::map<std::string,int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<std::endl;
//     }
//   return s;
// }

// std::ostream& operator<<(std::ostream &s, const std::map<std::string,double>& A)
// {
//   for(std::map<std::string,double>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<std::endl;
//     }
//   return s;
// }

// std::ostream& operator<<(std::ostream &s, const std::map<std::string,std::string>& A)
// {
//   for(std::map<std::string,std::string>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << "\t"<<p->first<<"\t --> "<<p->second<<std::endl;
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// std::ostream& operator<<(std::ostream &s, const std::set<std::string>& A)
// {
//   for(std::set<std::string>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << *p << " ";
//     }
//   return s;
// }

// std::ostream& operator<<(std::ostream &s, const std::set<int>& A)
// {
//   for(std::set<int>::const_iterator p=A.begin();p!=A.end();p++)
//     {
//       s << *p << " ";
//     }
//   return s;
// }

// /*-------------------------------------------------------------------*/

// std::ostream& operator<<(std::ostream& s, const std::vector<std::pair<int,int> >& A)
// {
//   typedef std::vector<std::pair<int,int> >::const_iterator it;
//   for(it p=A.begin();p!=A.end();p++)
//     {
//       s << "( "<<p->first <<" "<<p->second<<" )   ";
//     }
//   return s;
// }

void write_data(const CompVector<double>& v,std::ostream& s)
{
  s << v.n() << std::endl
    << v.ncomp() << std::endl
    << v << std::endl;
}

void read_data(CompVector<double>& v, std::istream& s)
{
  int n;
  int ncomp;
  s >> n;
  s >> ncomp;
  v.ncomp()=ncomp;
  v.resize(n,0);
  s >> v;
}

template<class T>
void write_data(const std::set<T>& v,std::ostream& s)
{
  s << v.size() << std::endl << v;
}

template<class T>
void read_data(std::set<T>& v, std::istream& s)
{
  size_t n;
  s >> n;
  for (int i=0;i<n;++i)
    {
      T a;
      s >> a;
      v.insert(a);
    }
}




void write_data(const int& v,std::ostream& s) { s << v; }
void read_data(int& v, std::istream& s)        { s >> v; }

template<int N,class T>
void write_data(const fixarray<N,T>& v,std::ostream& s) { s << v; }
template<int N,class T>
void read_data(fixarray<N,T>& v, std::istream& s)        { s >> v; }

void write_data(const double& v,std::ostream& s) { s << v; }
void read_data(double& v, std::istream& s)        { s >> v; }

void write_data(const std::string& v,std::ostream& s) { s << v; }
void read_data(std::string& v, std::istream& s)        { s >> v; }

template<class T>
void write_data(const std::vector<T>& v,std::ostream& s)
{
  s << v.size() << std::endl;
  for (int i=0;i<v.size();++i)
    {
      write_data(v[i],s);
      s << std::endl;
    }
  s << std::endl;
}


template<class T>
void read_data(std::vector<T>& v, std::istream& s)
{
  size_t n;
  s >> n;
  if(v.size()!=n) v.resize(n);
  for (int i=0;i<n;++i)
    read_data(v[i],s);
}

template void write_data(const fixarray<2,int>&,std::ostream& );
template void write_data(const fixarray<3,int>&,std::ostream& );
template void write_data(const fixarray<4,int>&,std::ostream& );
template void write_data(const fixarray<8,int>&,std::ostream& );
template void write_data(const fixarray<9,int>&,std::ostream& );

template void write_data(const std::vector<int>&,std::ostream& );
template void write_data(const std::vector<Vertex<2> >&,std::ostream& );
template void write_data(const std::vector<Vertex<3> >&,std::ostream& );
template void write_data(const std::vector<nvector<int> >&,std::ostream& );
template void write_data(const std::vector<fixarray<4,int> >&,std::ostream& );
template void write_data(const std::vector<fixarray<8,int> >&,std::ostream& );
template void write_data(const std::vector<std::set<int> >&,std::ostream& );

template void read_data(fixarray<2,int>&, std::istream& );
template void read_data(fixarray<3,int>&, std::istream& );
template void read_data(fixarray<4,int>&, std::istream& );
template void read_data(fixarray<8,int>&, std::istream& );
template void read_data(fixarray<9,int>&, std::istream& );

template void read_data(std::vector<int>&, std::istream& );
template void read_data(std::vector<Vertex<2> > &, std::istream& );
template void read_data(std::vector<Vertex<3> > &, std::istream& );
template void read_data(std::vector<nvector<int> >&, std::istream& );
template void read_data(std::vector<fixarray<4,int> >&, std::istream& );
template void read_data(std::vector<fixarray<8,int> >&, std::istream& );
template void read_data(std::vector<std::set<int> >&, std::istream& );
template void read_data(std::set<int> &, std::istream& );


