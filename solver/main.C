#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "eigen3/Eigen/Eigen"
#include "Sparsezzw.H"

using namespace std;
using namespace Eigen;

#include "define_function.H"

int main(int argc, char *argv[])
{
	cout<<"start !------------"<<endl;

    #include "createFields.H"

   if(SorT)
   {
   #include "SIMPLE_loop.H"
   }
   else 
   {
   #include "SIMPLE_t.H"	   
   }

   if(writet<0)
   {
   ofstream dataFileu;
   dataFileu.open("u.txt",ofstream::app);
        dataFileu<<u<<endl;
   dataFileu.close();

   ofstream dataFilev;
   dataFilev.open("v.txt",ofstream::app);
        dataFilev<<v<<endl;
   dataFilev.close();

   ofstream dataFilep;
   dataFilep.open("p.txt",ofstream::app);
        dataFilep<<p<<endl;
   dataFilep.close();
   }
    cout << "END !------------" << endl;

    return 0;
}
