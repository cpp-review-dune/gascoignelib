#ifndef __fixarray_h
#define __fixarray_h

#include  <stdlib.h> 
#include  <iostream> 
#include  <iterator> 
#include  <algorithm> 

/*-------------------------------------------------*/

template<int N, class T>
class fixarray
{

public:

  typedef  T*        iterator;
  typedef  const T*  const_iterator;
  
 protected:
  
  T  val[N];
  
  void array_copy(const_iterator q)
    {
      iterator       p(begin());
      const_iterator pe(end());
      while(p!=pe)  *p++ = *q++;
    }
  
 public:
  
  fixarray<N,T>()      { init(T());}
  fixarray<N,T>(const T& d) { init(d);}
  fixarray<N,T>(const fixarray<N,T>& v)
    {
      init(T());
      array_copy(v.begin());
      //copy(v.begin(),v.end(),begin());
    }
  fixarray(const_iterator b)
    {
      init(T());
      array_copy(b);
    }
  
  virtual ~fixarray()
    {
//       Destroy(begin(),end());
    }
  
  void init(const T& d)
    {
      // Braucht man das wirklich ???
//       for(int i=0;i<N;i++)  construct(&(val[i]),d);
      for(int i=0;i<N;i++)  val[i]=d;
    }
  
  const T*  begin() const { return &(val[0]);}
  const T*  end  () const { return &(val[0])+N;}
  T*        begin()       { return &(val[0]);}
  T*        end  ()       { return &(val[0])+N;}
  
  size_t   size()            const { return N;}
  const T& operator[](int i) const { return val[i];}
  T&       operator[](int i)       { return val[i];}
  
  fixarray<N,T>& operator=(const T& d) 
    {
      iterator  p(end());
      while(p>begin()) *--p = d;
      return *this;
    } 
  
  fixarray<N,T>& operator=(const fixarray<N,T>& v) 
    {
      iterator        p(begin());
      const_iterator  q(v.begin());
      while(p<end()) *p++ = *q++;
      return *this;
    } 
  
  bool operator<(const fixarray<N,T>& v) const
    {
      const_iterator  p(  begin());
      const_iterator  q(v.begin());
      while(p<end())
	{
	  if (*p<*q) return 1;
	  if (*q<*p) return 0;
	  p++; q++;
	}
      return 0;
    }
  bool operator!=(const fixarray<N,T>& v) const
    {
      const_iterator  p(  begin());
      const_iterator  q(v.begin());
      while(p<end())
	{
	  if (*p!=*q) return 1;
	  p++; q++;
	}
      return 0;
    }
  
    
  std::ostream& put(std::ostream &s) const
    {
      copy(begin(),end(),std::ostream_iterator<T>(s," "));
      return s;
    }  
  
  std::istream& get(std::istream &s)
    {
      for(fixarray<N,T>::iterator p = begin();p!=end();p++) s >> *p;
      return s;
    }
  
  void read_data(std::istream& s)
    {
      size_t n;
      s >> n;
      if(size()!=n) 
	{
	  std::cerr << "read_data(): wrong size in fixarray" << N << " " << n << std::endl;
	  exit(1);
	}
      s >> *this;
    }
  
  void write_data(std::ostream& s) const
    {
      s << size() << std::endl;
      s << *this;
    }
};

/*-------------------------------------------------*/

template<int N, class T>
bool operator==(const fixarray<N,T>& x, const fixarray<N,T>& y) 
{
  return std::equal(x.begin(), x.end(), y.begin());
}

/*-------------------------------------------------*/

class fixarrayHash
{
 public:
  template<int N, class T>
    int operator()(const fixarray<N,T>& h) const { return (int) h[0];}
};


template<int N,class T>
std::ostream& operator<<(std::ostream &s, const fixarray<N,T>& A);
template<int N,class T>
std::istream& operator>>(std::istream &s, fixarray<N,T>& A);



#endif


