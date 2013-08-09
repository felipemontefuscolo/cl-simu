#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>

double pi=3.14;

using namespace std;

bool is_operator(char c)
{
  if(c=='*' || c=='/' || c=='+' || c=='-')
    return true;
  else
    return false;
};

int is_par(char c)
{
  if(c=='(')
    return 1;
  if(c==')')
    return -1;
  return 0;
}

bool is_delim(char c)
{
  if(c=='[')
    return 1;
  if(c==']')
    return -1;
  return 0;
}

template<class T>
void findandreplace( T& source, const T& find, const T& replace )
{
size_t j;
  for (;(j = source.find( find )) != T::npos;)
  {
    source.replace( j, find.length(), replace );
  }
};

string next_term(string const&source)
{
  char c;
  int  np(0);   // num de parenteses 
  string term;
  bool open_p(false);
  bool is_op(false);
  
  for(int i=0;i!=source.size();i++)
  {
    c = source[i];          //cout << "c= "<<c << endl;
    np   += is_par(c);      //cout << "np= " <<np << endl;
    is_op = is_operator(c); //cout << "is_op= "<< is_op << endl;
    if (np>0)
      is_op = false;
    
    if (is_op || np < 0)
      break;
    
    term.append(1,c);
  }

  return term;
};

string prev_term(string const&source)
{
  char c;
  int  np(0);   // num de parenteses 
  string term;
  bool open_p(false);
  bool is_op(false);
  
  for(int i=source.size()-1;i>=0;i--)
  {
    c = source[i];          //cout << "c= "<<c << endl;
    np   -= is_par(c);      //cout << "np= " <<np << endl;
    is_op = is_operator(c); //cout << "is_op= "<< is_op << endl;
    if (np>0)
      is_op = false;
    
    if (is_op || np < 0)
      break;
    
    term.append(1,c);
  }

  reverse(term.begin(), term.end());
  

  return term;
};

void fix_pot(string &exp)
{
  char    c;
  size_t  t, j;  
  
  for (;(j = exp.find( "^" )) != string::npos;)
  {
    //j = exp.find( "^" );
    string lhs, rhs, l, r;
    rhs = next_term(  exp.substr(j+1,exp.size()-j)  );
    lhs = prev_term(  exp.substr(0, j)              );

    r = rhs;
    l = lhs;
    if (r[0]=='(') {r.erase(r.begin());r.erase(r.end()-1);};
    if (l[0]=='(') {l.erase(l.begin());l.erase(l.end()-1);};
  
    if (lhs==string("%e"))
      findandreplace(exp, lhs+"^"+rhs, "exp("+r+")");
    else
      findandreplace(exp, lhs+"^"+rhs, "pow("+l+","+r+")");
  }
  
}

vector<string> find_exp(string argument)
{
  vector<string> temp;
  vector<string> exps;
  size_t         j,k=0;
  
  // detectando o que está entre colchetes  
  j = argument.find( '[' );
  if (j==string::npos) {
    temp.push_back(argument);
  }
  for (;(j = argument.find( '[',k )) != string::npos;)
  {
    k = argument.find( ']',k+1 );
    temp.push_back( argument.substr(j+1,k-j-1) );
  }
  
  // detectando o que está separado por vírgulas
  j = argument.find( ',' );
  if (j==string::npos) {
    return temp;
  }
  for (int i = 0; i < temp.size(); i++)
  {
    k=0;
    for (;(j = temp[i].find( ',',k )) != string::npos;)
    {
      exps.push_back( temp[i].substr(k,j-k) );
      k = j + 1;
    }
    exps.push_back( temp[i].substr(k,temp[i].size()-k) );
  }
  
  return exps;
}

int main(int argc, char *argv[])
{
  string         argument;
  vector<string> exps;
  
  if (argc < 2)
  {
    cout << "uso:" << endl;
    cout << argv[0] << " \"expressao\"" << endl;
    return 0;
  }
    
  argument.assign(argv[1]);
  
  // elimina espaços em branco
  findandreplace(argument, string(" "), string(""));
  
  exps = find_exp(argument);
  
  cout << endl;
  for (int i = 0; i < exps.size(); i++)
  {
    fix_pot(exps[i]);
    cout << exps[i] << ";" << endl;
  }
  cout << endl;

  return 0;
};


//\[\left( -48\,{x}^{2}+48\,x-8\right) \,{y}^{3}+\left( 72\,{x}^{2}-72\,x+12\right) \,{y}^{2}+\left( -24\,{x}^{4}+48\,{x}^{3}-48\,{x}^{2}+24\,x-4\right) \,y+12\,{x}^{4}-24\,{x}^{3}+12\,{x}^{2}-2\,x+1\]






















  K[1][i]  : at(H[i], [l1=0  ,  l2=1,   l3=0  ]),
  K[2][i]  : at(H[i], [l1=1/3,  l2=2/3, l3=0  ]),
  K[3][i]  : at(H[i], [l1=2/3,  l2=1/3, l3=0  ]),
  K[4][i]  : at(H[i], [l1=1  ,  l2=0,   l3=0  ]),
  K[5][i]  : at(H[i], [l1=0  ,  l2=2/3, l3=1/3]),
  K[6][i]  : at(H[i], [l1=1/3,  l2=1/3, l3=1/3]),
  K[7][i]  : at(H[i], [l1=2/3,  l2=0  , l3=1/3]),
  K[8][i]  : at(H[i], [l1=0  ,  l2=1/3, l3=2/3]),
  K[9][i]  : at(H[i], [l1=1/3,  l2=0  , l3=2/3]),
  K[10][i] : at(H[i], [l1=0  ,  l2=0,   l3=1  ]),






  R[1][1]  : at((l1+l2)^2 * (l1-l2), [l1=0  ,  l2=1,   l3=0  ]),
  R[2][1]  : at((l1+l2)^2 * (l1-l2), [l1=1/3,  l2=2/3, l3=0  ]),
  R[3][1]  : at((l1+l2)^2 * (l1-l2), [l1=2/3,  l2=1/3, l3=0  ]),
  R[4][1]  : at((l1+l2)^2 * (l1-l2), [l1=1  ,  l2=0,   l3=0  ]),
  R[5][1]  : at((l1+l2)^2 * (l1-l2), [l1=0  ,  l2=2/3, l3=1/3]),
  R[6][1]  : at((l1+l2)^2 * (l1-l2), [l1=1/3,  l2=1/3, l3=1/3]),
  R[7][1]  : at((l1+l2)^2 * (l1-l2), [l1=2/3,  l2=0  , l3=1/3]),
  R[8][1]  : at((l1+l2)^2 * (l1-l2), [l1=0  ,  l2=1/3, l3=2/3]),
  R[9][1]  : at((l1+l2)^2 * (l1-l2), [l1=1/3,  l2=0  , l3=2/3]),
  R[10][1] : at((l1+l2)^2 * (l1-l2), [l1=0  ,  l2=0,   l3=1  ])







