double A(double f,double d)
{
    return 1;
    //Can modify to set inviscid flow, for the expression of infinite Reynolds number
}

double psi(double r)
{
    double mini=std::min(2.0*r,(1.0+3.0*r)/4.0);
    mini = std::min(mini, (3.0 + r) / 4.0);
    mini = std::min(mini, 2.0);
    double psi=std::max(0.0,mini);
    //Van Leer, Adam ~~
    //double psi=std::max(  0.0 , ( r+pow(r,2) ) / ( 1.0+pow(r,2) )  );
    return psi;
}


