#ifndef SPARSEZZW_H
#define SPARSEZZW_H

//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <string>
//#include <sstream>
#include "eigen3/Eigen/Eigen"

class Sparsezzw
{
//private:
    // Private Member Functions   
public:
	Eigen::MatrixXd N5;
	int Ni;
    // Constructor
        Sparsezzw(const int ni, const int nj);
        
    // Member functions
        void input(const int i,const int j,const double iinput);
        double output(const int i,const int j) const;

    //- Destructor
        virtual ~Sparsezzw();
};


#endif
