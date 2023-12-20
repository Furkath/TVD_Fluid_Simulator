#include <iostream>
#include "Sparsezzw.H"
//#include <cmath>
//#include <algorithm>
//#include <time.h>
using namespace std;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Sparsezzw::Sparsezzw(const int ni, const int nj)
{
    N5 = Eigen::MatrixXd::Zero(ni*nj,5);
    Ni=ni;
}

// Public Member Functions
void Sparsezzw::input(const int i, const int j, const double iinput)
{
    //cout<<i<<"  sparse i  "<<j<<endl;
    if (i==j)
    {
	//    cout<<"is i  "<<j<<endl;
        N5(i, 2) = iinput;
    }
    else if (j==i-1)
    {
	//    cout<<"is i  "<<j<<endl;
        N5(i, 1) = iinput;
    }
    else if (j==i+1)
    {
	//    cout<<"is i  "<<j<<endl;
        N5(i, 3) = iinput;
    }
    else if (j==i-Ni)
    {
	//    cout<<"is i  "<<j<<endl;
        N5(i, 0) = iinput;
    }
    else if (j==i+Ni)
    {
	//    cout<<"this ?  "<<j<<endl;
        N5(i, 4) = iinput;
        //cout<<"~<>~"<<endl;
    }
    else
    {
        cerr<<"j index not considered!"<<endl;
    }
    //cout<<"into "<<endl;
}
double Sparsezzw::output(const int i, const int j) const
{
    //cout<<"into "<<endl;
    //cout << i << "  sparse O  " << j << endl;
    if (i == j)
    {
	//    cout<<"is o  "<<j<<endl;
        return N5(i, 2);
    }
    else if (j == i - 1)
    {
	//    cout<<"is o  "<<j<<endl;
        return N5(i, 1);
    }
    else if (j == i + 1)
    {
	//    cout<<"is o  "<<j<<endl;
        return N5(i, 3);
    }
    else if (j == i - Ni)
    {
	//    cout<<"is o  "<<j<<endl;
        return N5(i, 0);
    }
    else if (j == i + Ni)
    {
	//    cout<<"is o  "<<j<<endl;
        return N5(i, 4);
    }
    else
    {
        cerr << "j index not considered!" << endl;
    }
}

Sparsezzw::~Sparsezzw()
{
}
