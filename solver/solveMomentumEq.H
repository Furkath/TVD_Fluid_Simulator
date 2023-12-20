{
    double De = mu * dy / dx;
    double Dw = mu * dy / dx;
    double Dn = mu * dx / dy;
    double Ds = mu * dx / dy;
    double Fe;
    double Fw;
    double Fn;
    double Fs;
    double aE;
    double aW;
    double aN;
    double aS;
    double aP;
    double pressure_term;


    u_star = MatrixXd::Zero(ii + 1, jj);
    d_u = MatrixXd::Zero(ii + 1, jj);

    for (int i=1;i<ii;i++)
    {
        for (int j=1;j<jj-1;j++)
        {
            Fe = 0.5 * rho * dy * (u(i + 1, j) + u(i, j));
            Fw = 0.5 * rho * dy * (u(i - 1, j) + u(i, j));
            Fn = 0.5 * rho * dx * (v(i, j + 1) + v(i - 1, j + 1));
            Fs = 0.5 * rho * dx * (v(i, j) + v(i - 1, j));

            aE = De * A(Fe, De) + max(-Fe, 0.0);
            aW = Dw * A(Fw, Dw) + max(Fw, 0.0);
            aN = Dn * A(Fn, Dn) + max(-Fn, 0.0);
            aS = Ds * A(Fs, Ds) + max(Fs, 0.0);
            aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);

            pressure_term = (p(i - 1, j) - p(i, j)) * dy;

            u_star(i, j) = alphaS / aP * ((aE * u(i + 1, j) + aW * u(i - 1, j) + aN * u(i, j + 1) + aS * u(i, j - 1)) + pressure_term) + (1 - alphaS) * u(i, j);

            d_u(i, j) = alphaS * dy / aP;
        }
    }
    int j=0;
    for(int i=1;i<ii;i++)
    {
        Fe = 0.5 * rho * dy * (u(i + 1, j) + u(i, j));
        Fw = 0.5 * rho * dy * (u(i - 1, j) + u(i, j));
        Fn = 0.5 * rho * dx * (v(i, j + 1) + v(i - 1, j + 1));
        Fs = 0;

        aE = De * A(Fe, De) + max(-Fe, 0.0);
        aW = Dw * A(Fw, Dw) + max(Fw, 0.0);
        aN = Dn * A(Fn, Dn) + max(-Fn, 0.0);
        aS = 0;
        aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);
        d_u(i, j) = alphaS * dy / aP;
    }
    j=jj-1;
    for (int i = 1; i < ii; i++)
    {
        Fe = 0.5 * rho * dy * (u(i + 1, j) + u(i, j));
        Fw = 0.5 * rho * dy * (u(i - 1, j) + u(i, j));
        Fn = 0;
        Fs = 0.5 * rho * dx * (v(i, j) + v(i - 1, j));

        aE = De * A(Fe, De) + max(-Fe, 0.0);
        aW = Dw * A(Fw, Dw) + max(Fw, 0.0);
        aN = 0;
        aS = Ds * A(Fs, Ds) + max(Fs, 0.0);
        aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);
        d_u(i, j) = alphaS * dy / aP;
    }
    for (int J=0;J<jj;J++)
    {
        u_star(0,J)=-u_star(1,J);
        u_star(ii,J)=-u_star(ii-1,J);
    }
    for (int I=0;I<ii+1;I++)
    {
        u_star(I,0)=0.0;
        u_star(I,jj-1)=v_top;
    }

    v_star = MatrixXd::Zero(ii, jj+1);
    d_v = MatrixXd::Zero(ii, jj+1);
    De = mu * dy / dx;
    Dw = mu * dy / dx;
    Dn = mu * dx / dy;
    Ds = mu * dx / dy;
    for(int i=1;i<ii-1;i++)
    {
        for (int j=1;j<jj;j++)
        {
            Fe = 0.5 * rho * dy * (u(i + 1, j) + u(i + 1, j - 1));
            Fw = 0.5 * rho * dy * (u(i, j) + u(i, j - 1));
            Fn = 0.5 * rho * dx * (v(i, j) + v(i, j + 1));
            Fs = 0.5 * rho * dx * (v(i, j - 1) + v(i, j));

            aE = De * A(Fe, De) + max(-Fe, 0.0);
            aW = Dw * A(Fw, Dw) + max(Fw, 0.0);
            aN = Dn * A(Fn, Dn) + max(-Fn, 0.0);
            aS = Ds * A(Fs, Ds) + max(Fs, 0.0);
            aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);

            pressure_term = (p(i, j - 1) - p(i, j)) * dx;

            v_star(i, j) = alphaS / aP * ((aE * v(i + 1, j) + aW * v(i - 1, j) + aN * v(i, j + 1) + aS * v(i, j - 1)) + pressure_term) + (1 - alphaS) * v(i, j);

            d_v(i, j) = alphaS * dx / aP;
        }
    }
    int i=0;
    for (int j=1;j<jj;j++)
    {
        Fe = 0.5 * rho * dy * (u(i + 1, j) + u(i + 1, j - 1));
        Fw = 0;
        Fn = 0.5 * rho * dx * (v(i, j) + v(i, j + 1));
        Fs = 0.5 * rho * dx * (v(i, j - 1) + v(i, j));

        aE = De * A(Fe, De) + max(-Fe, 0.0);
        aW = 0;
        aN = Dn * A(Fn, Dn) + max(-Fn, 0.0);
        aS = Ds * A(Fs, Ds) + max(Fs, 0.0);
        aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);
        d_v(i, j) = alphaS * dx / aP;
    }
    i=ii-1;
    for (int j = 1; j < jj; j++)
    {
        Fe = 0;
        Fw = 0.5 * rho * dy * (u(i, j) + u(i, j - 1));
        Fn = 0.5 * rho * dx * (v(i, j) + v(i, j + 1));
        Fs = 0.5 * rho * dx * (v(i, j - 1) + v(i, j));

        aE = 0;
        aW = Dw * A(Fw, Dw) + max(Fw, 0.0);
        aN = Dn * A(Fn, Dn) + max(-Fn, 0.0);
        aS = Ds * A(Fs, Ds) + max(Fs, 0.0);
        aP = aE + aW + aN + aS + (Fe - Fw) + (Fn - Fs);
        d_v(i, j) = alphaS * dx / aP;
    }
    for (int J=0;J<jj+1;J++)
    {
        v_star(0,J)=0.0;
        v_star(ii-1,J)=0.0;
    }
    for (int I=0;I<ii;I++)
    {
        v_star(I,0)=-v_star(I,1);
        v_star(I,jj)=-v_star(I,jj-1);
    }
}


































/*
//Solve U E.q. with Gauss-Siedel iteration
n_U_loop=0;
while (true)
{
    n_U_loop++;
    double old= 0;
    double residual=0;
    UorV=1;
    #include "updateCoefficients.H"
    for (int i=0;i<N-1;i++)
    {
        for (int j=0;j<N;j++)
        {
            old=U(i,j);
            U(i, j) = U(i,j)*(1-as) + as*1.0 / AUiJ(i, j) * (
                AUi_1J(i, j) *   (  i==0   ? 0 : U(i - 1, j)  )
                + AUi1J(i, j) *  (  i==N-2 ? 0 : U(i + 1, j)  )
                + AUiJ_1(i, j) * (  j==0   ? 0 : U(i, j - 1)  )
                + AUiJ1(i, j) *  (  j==N-1 ? 0 : U(i, j + 1)  )
                + SUiJ(i, j)
            );
	    //if(i==N-2) cout<<"Utop:  "<<U(N-2,j)<<endl;
            residual=max(residual,abs(old-U(i,j)));
        }
    }
    if (residual<tol) 
    {
	    cout<<"No. U loop:  "<<n_U_loop<<endl;
	    
	    break;
    }
}
// Solve V E.q. with Gauss-Siedel iteration
n_V_loop=0;
while (true)
{
    n_V_loop++;
    double old=0;
    double residual=0;
    UorV=0;
    #include "updateCoefficients.H"
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N-1;j++)
        {
            old=V(i,j);
            V(i, j) = V(i,j)*(1-as) + as*1.0 / AVIj(i, j)*(
                 AVI_1j(i, j) *   (  i==0   ? 0 : V(i - 1, j)  )
                 + AVI1j(i, j) *  (  i==N-1 ? 0 : V(i + 1, j)  )
                 + AVIj_1(i, j) * (  j==0   ? 0 : V(i, j - 1)  )
                 + AVIj1(i, j) *  (  j==N-2 ? 0 : V(i, j + 1)  )
                 + SVIj(i, j)
            );
            residual=max(residual,abs(old-V(i,j)));
        }
    }
    if (residual < tol)
    {
            cout<<"No. V loop:  "<<n_V_loop<<endl;
            break;
    }
}
*/