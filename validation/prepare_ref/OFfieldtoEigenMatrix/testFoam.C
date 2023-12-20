/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"
#include "eigen3/Eigen/Eigen"
#include <cmath>
using namespace Eigen;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   label cellSize = mesh.cells().size();
   double cellsize = cellSize;
   double nCells = std::sqrt(cellsize);

   int length=nCells;

   scalar x;
   scalar y;
   scalar ii;
   scalar jj;
   double i_;
   double j_;
   int i;
   int j;
   MatrixXd u(length,length);
   MatrixXd v(length,length);
   //Info<<length<<endl;
   forAll(U,celli)
   {
     x=mesh.C()[celli].x();
     y=mesh.C()[celli].y();
     //Info<<x<<"  "<<y<<endl;
     ii= (x-0.5/nCells)*nCells;
     jj= (y-0.5/nCells)*nCells;
     i_=ii;
     j_=jj;
     i=round(i_);
     j=round(j_);
     //Info<<ii<<" "<<jj<<endl;
     //Info<<i<<"  "<<j<<endl;
     u(i,j)=U[celli].x();
     v(i,j)=U[celli].y();
     //Info<<u(i,j)<<"  "<<v(i,j)<<endl;
   }
   //Info<<nCells<<endl;
   //cout<<u<<endl;
   //cout<<v<<endl;
   ofstream dataFilep;
   dataFilep.open("u_ref.txt",ofstream::app);
   //fstream file("result.txt", ios::out);
   //for (int i=0;i<0.5*u.size();i++)
   //{
        //dataFile<<std::log(k[i])<<" "<<std::log(E[i])<<endl;
        dataFilep<<u<<endl;
   //}
   dataFilep.close();

   //ofstream dataFilep;
   dataFilep.open("v_ref.txt",ofstream::app);
   //fstream file("result.txt", ios::out);
   //for (int i=0;i<0.5*u.size();i++)
   //{
        //dataFile<<std::log(k[i])<<" "<<std::log(E[i])<<endl;
        dataFilep<<v<<endl;
   //}
   dataFilep.close();


    return 0;
}


// ************************************************************************* //
