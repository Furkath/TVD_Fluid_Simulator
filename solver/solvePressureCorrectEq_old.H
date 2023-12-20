{
    VectorXd bp = Eigen::VectorXd::Zero(ii * jj);
    int dist=ii;
    int location;
    double aE;
    double aW;
    double aN;
    double aS;
    double aP;
    for (int i = 0; i < ii; i++)
    {
        for (int j=0;j<jj;j++)
        {
            location=i+(j-0)*dist;
            bp(location) = rho * (u_star(i, j) * dy - u_star(i + 1, j) * dy + v_star(i, j) * dx - v_star(i, j + 1) * dx);
        }
    }
    bp(0)=0.0;
    //cout<<bp<<endl;
    Sparsezzw Ap(ii,jj);
    //cout<<"Here"<<endl;
    for (int j=0;j<jj;j++)
    {
        for (int i=0;i<ii;i++)
        {
	    //	cout<<i<<"  "<<j<<endl;
            location = i + (j - 0) * dist;
            aE = 0;
            aW = 0;
            aN = 0;
            aS = 0;
            if(i==0 && j==0)
            {
                Ap.input(location,location,1);
            }
            else if(i==ii-1 && j==0)
            {
		//    cout<<"AAAA"<<endl;
                Ap.input(location, location - 1, -rho * d_u(i, j) * dy);
                aW = -Ap.output(location, location - 1);
		//cout<<"The 1"<<endl;
		//cout<<location<<" + "<<dist<<"ii jj"<<ii<<"  "<<jj<<endl;
		//cout<<"!!!!!!!!!!!!!!!"<<endl;
		//cout<<location<<endl;
		//cout<<location+dist<<endl;
		//cout<<Ap.N5<<endl;

                Ap.input(location, location+dist ,  -rho * d_v(i, j + 1) * dx);
		//cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                aN = -Ap.output(location, location + dist);
		//cout<<"The 2"<<endl;
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
		//cout<<"BBBBBB"<<endl;
            }
            else if(i==0 && j==jj-1)
            {
                Ap.input(location, location + 1, -rho * d_u(i + 1, j) * dy);
                aE = -Ap.output(location, location + 1);
                Ap.input(location, location - dist, -rho * d_v(i, j) * dx);
                aS = -Ap.output(location, location - dist);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            } 
            else if(i==ii-1 && j==jj-1)
            {
                Ap.input(location, location - 1, -rho * d_u(i, j) * dy);
                aW = -Ap.output(location, location - 1);
                Ap.input(location, location - dist, -rho * d_v(i, j) * dx);
                aS = -Ap.output(location, location - dist);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            }
            else if(i==0)
            {
                Ap.input(location, location + 1, -rho * d_u(i + 1, j) * dy);
                aE = -Ap.output(location, location + 1);
                Ap.input(location, location + dist, -rho * d_v(i, j + 1) * dx);
                aN = -Ap.output(location, location + dist);
                Ap.input(location, location - dist, -rho * d_v(i, j) * dx);
                aS = -Ap.output(location, location - dist);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            }
            else if(j==0)
            {
                Ap.input(location, location + 1, -rho * d_u(i + 1, j) * dy);
                aE = -Ap.output(location, location + 1);
                Ap.input(location, location + dist, -rho * d_v(i, j + 1) * dx);
                aN = -Ap.output(location, location + dist);
                Ap.input(location, location - 1, -rho * d_u(i, j) * dy);
                aW = -Ap.output(location, location - 1);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            }
            else if(i==ii-1)
            {
                Ap.input(location, location + dist, -rho * d_v(i, j + 1) * dx);
                aN = -Ap.output(location, location + dist);
                Ap.input(location, location - dist, -rho * d_v(i, j) * dx);
                aS = -Ap.output(location, location - dist);
                Ap.input(location, location - 1, -rho * d_u(i, j) * dy);
                aW = -Ap.output(location, location - 1);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            }
            else if(j==jj-1)
            {
                Ap.input(location, location + 1, -rho * d_u(i + 1, j) * dy);
                aE = -Ap.output(location, location + 1);
                Ap.input(location, location - dist, -rho * d_v(i, j) * dx);
                aS = -Ap.output(location, location - dist);
                Ap.input(location, location - 1, -rho * d_u(i, j) * dy);
                aW = -Ap.output(location, location - 1);
                aP = aE + aN + aW + aS;
                Ap.input(location, location, aP);
            }
            else
            {
                Ap.input(location,location-1,-rho*d_u(i,j)*dy);  
                aW = -Ap.output(location,location-1);

                Ap.input(location,location+1,-rho*d_u(i+1,j)*dy);     
                aE = -Ap.output(location,location+1);

                Ap.input(location,location-dist,-rho*d_v(i,j)*dx);  
                aS = -Ap.output(location,location-dist);

                Ap.input(location,location+dist,-rho*d_v(i,j+1)*dx);     
                aN = -Ap.output(location,location+dist);
                
                aP = aE + aN + aW + aS;
                Ap.input(location,location,aP);
            }
        }
    }
    //cout<<Ap.N5<<endl;
p=p_star;
p_prime = Eigen::MatrixXd::Zero(ii, jj);

int N=ii*jj;
VectorXd x = Eigen::VectorXd::Zero(N);
//cout<<"Here 2"<<endl;


VectorXd d = Ap.N5.col(2);
VectorXd e = Ap.N5.col(3).head(N-1);
VectorXd f = Eigen::VectorXd::Zero(N-2);//Ap.N5.col(4).head(ii*jj-2);
VectorXd h = Ap.N5.col(1);
h(0)=0.0;
VectorXd g = Eigen::VectorXd::Zero(N);

VectorXd alpha = Eigen::VectorXd::Zero(N);
VectorXd gamma = Eigen::VectorXd::Zero(N-1);
VectorXd delta = Eigen::VectorXd::Zero(N-2);
VectorXd beta = Eigen::VectorXd::Zero(N);
VectorXd c = Eigen::VectorXd::Zero(N);
VectorXd z = Eigen::VectorXd::Zero(N);

alpha(0) = d(0);
gamma(0) = e(0) / alpha(0);
delta(0) = f(0) / alpha(0);

beta(1)=h(1);
alpha(1) = d(1) - beta(1) * gamma(0);
gamma(1) = (e(1) - beta(1) * delta(0)) / alpha(1);
delta(1) = f(1) / alpha(1);
for (int k=2;k<N-2;k++)
{
    beta(k) = h(k) - g(k) * gamma(k - 2);
    alpha(k) = d(k) - g(k) * delta(k - 2) - beta(k) * gamma(k - 1);
    gamma(k) = (e(k) - beta(k) * delta(k - 1)) / alpha(k);
    delta(k) = f(k) / alpha(k);
}
beta(N - 2) = h(N - 2) - g(N - 2) * gamma(N - 4);
alpha(N - 2) = d(N - 2) - g(N - 2) * delta(N - 4) - beta(N - 2) * gamma(N - 3);
gamma(N - 2) = (e(N - 2) - beta(N - 2) * delta(N - 3)) / alpha(N - 2);
beta(N-1) = h(N-1) - g(N-1) * gamma(N - 3);
alpha(N-1) = d(N-1) - g(N-1) * delta(N - 3) - beta(N-1) * gamma(N - 2);

c(0) = bp(0) / alpha(0);
c(1) = (bp(1) - beta(1) * c(0)) / alpha(1);
for (int k=2;k<N;k++)
{
    c(k) = (bp(k) - g(k) * c(k - 2) - beta(k) * c(k - 1)) / alpha(k);
}
//c=z.cwiseQuotient(alpha);
x(N-1)=c(N-1);
x(N-2)=c(N-2)-gamma(N-2)*x(N-1);
for (int k=N-3;k>=0;k--)
{
    x(k) = c(k) - gamma(k) * x(k + 1) - delta(k) * x(k + 2);
}


int pos=0;
for (int j=0;j<jj;j++)
{
    for (int i=0;i<ii;i++)
    {
        p_prime(i, j) = x(pos);
        pos++;
        p(i, j) = p_star(i, j) + alphaP * p_prime(i, j);
    }
}
p(0,0)=0.0;
}





/*
Pc = MatrixXd::Zero(N - 1, N - 1);
n_Pc_loop=0;
while(true)
{ 
n_Pc_loop++;
double old = 0;
double residual = 0;
//Consider B.C. as zero-gradient
        /*
        old=Pc(0,0);
        Pc(0,0) = Pc(0,0)*(1-as)+
        as*1.0/(  /*1.0/AUiJ(0,0)+*///1.0/AUiJ(0,0+1)/*+1.0/AVIj(0,0)*/+1.0/AVIj(0+1,0)  )*( 
        // U(0,0)   +1.0/AUiJ(0,0)*   Pc(0,0-1) 
        /* -U(0,0+1) +1.0/AUiJ(0,0+1)* Pc(0,0+1)  
        //+V(0,0)   +1.0/AVIj(0,0)*   Pc(0-1,0)
        -V(0+1,0) +1.0/AVIj(0+1,0)* Pc(0+1,0) 
        );
        */
/*
        //To fix a reference value, otherwise the pc soulution would be undetermined,
        //as adding a const for all points would still satisfy the equation.
        //Physically, this is like to treat point(0,0)'s pressure as a refernece if the
        //cavaty is isolated.
        //Mathematically, this is to satisfy the Scarborough Criterion to make the coefficient
        //matrix bounded and diagonable, thus the iteration can converge.
        //And more weight of the central coefficient, faster the converge speed,
        //so here directly set the fist equation as the value assignment form.
        Pc(0,0)=0.0;
        //residual=max(residual,abs(Pc(0,0)-old));
        
        old=Pc(N-2,N-2);
        Pc(N-2,N-2) = Pc(N-2,N-2)*(1-as)+
        as*1.0/(  1.0/AUiJ(N-2,N-2)/*+1.0/AUiJ(N-2,N-2+1)*//*+1.0/AVIj(N-2,N-2)/*+1.0/AVIj(N-2+1,N-2)*/  /*)*(
         U(N-2,N-2)   +1.0/AUiJ(N-2,N-2)*   Pc(N-2,N-2-1)
        //-U(N-2,N-2+1) +1.0/AUiJ(N-2,N-2+1)* Pc(N-2,N-2+1) 
        +V(N-2,N-2)   +1.0/AVIj(N-2,N-2)*   Pc(N-2-1,N-2)
        //-V(N-2+1,N-2) +1.0/AVIj(N-2+1,N-2)* Pc(N-2+1,N-2) 
        );
        residual=max(residual,abs(Pc(N-2,N-2)-old));

        old=Pc(0,N-2);
        Pc(0,N-2) = Pc(0,N-2)*(1-as)+
        as*1.0/(  1.0/AUiJ(0,N-2)/*+1.0/AUiJ(0,N-2+1)+1.0/AVIj(0,N-2)*//*+1.0/AVIj(0+1,N-2)  )*(
         U(0,N-2)   +1.0/AUiJ(0,N-2)*   Pc(0,N-2-1)
        //-U(0,N-2+1) +1.0/AUiJ(0,N-2+1)* Pc(0,N-2+1) 
        //+V(0,N-2)   +1.0/AVIj(0,N-2)*   Pc(0-1,N-2)
        -V(0+1,N-2) +1.0/AVIj(0+1,N-2)* Pc(0+1,N-2) 
        );
        residual=max(residual,abs(Pc(0,N-2)-old));

        old=Pc(N-2,0);
        Pc(N-2,0) = Pc(N-2,0)*(1-as)+
        as*1.0/(  /*1.0/AUiJ(N-2,0)+*//*1.0/AUiJ(N-2,0+1)+1.0/AVIj(N-2,0)/*+1.0/AVIj(N-2+1,0)*/ /* )*(
        // U(N-2,0)   +1.0/AUiJ(N-2,0)*   Pc(N-2,0-1)
        -U(N-2,0+1) +1.0/AUiJ(N-2,0+1)* Pc(N-2,0+1) 
        +V(N-2,0)   +1.0/AVIj(N-2,0)*   Pc(N-2-1,0)
        //-V(N-2+1,0) +1.0/AVIj(N-2+1,0)* Pc(N-2+1,0) 
        );
        residual=max(residual,abs(Pc(N-2,0)-old));
//continue B.C. as zero-gradient
for (int j=1;j<N-2;j++)
{
        old=Pc(0,j);
        Pc(0,j) = Pc(0,j)*(1-as)+
        as*1.0/(  1.0/AUiJ(0,j)+1.0/AUiJ(0,j+1)/*+1.0/AVIj(0,j)*//*+1.0/AVIj(0+1,j)  )*(
         U(0,j)   +1.0/AUiJ(0,j)*   Pc(0,j-1)
        -U(0,j+1) +1.0/AUiJ(0,j+1)* Pc(0,j+1) 
        //+V(0,j)   +1.0/AVIj(0,j)*   Pc(0-1,j)
        -V(0+1,j) +1.0/AVIj(0+1,j)* Pc(0+1,j) 
        );
        residual=max(residual,abs(Pc(0,j)-old));

        old=Pc(N-2,j);
        Pc(N-2,j) = Pc(N-2,j)*(1-as)+
        as*1.0/(  1.0/AUiJ(N-2,j)+1.0/AUiJ(N-2,j+1)+1.0/AVIj(N-2,j)/*+1.0/AVIj(N-2+1,j)*/ /* )*(
         U(N-2,j)   +1.0/AUiJ(N-2,j)*   Pc(N-2,j-1)
        -U(N-2,j+1) +1.0/AUiJ(N-2,j+1)* Pc(N-2,j+1) 
        +V(N-2,j)   +1.0/AVIj(N-2,j)*   Pc(N-2-1,j)
        //-V(N-2+1,j) +1.0/AVIj(N-2+1,j)* Pc(N-2+1,j) 
        );
        residual=max(residual,abs(Pc(N-2,j)-old));
}
for (int i=1;i<N-2;i++)
{
        old=Pc(i,0);
        Pc(i,0) = Pc(i,0)*(1-as)+
        as*1.0/(  /*1.0/AUiJ(i,0)+*//*1.0/AUiJ(i,0+1)+1.0/AVIj(i,0)+1.0/AVIj(i+1,0)  )*(
        // U(i,0)   +1.0/AUiJ(i,0)*   Pc(i,0-1)
        -U(i,0+1) +1.0/AUiJ(i,0+1)* Pc(i,0+1) 
        +V(i,0)   +1.0/AVIj(i,0)*   Pc(i-1,0)
        -V(i+1,0) +1.0/AVIj(i+1,0)* Pc(i+1,0) 
        );
        residual=max(residual,abs(Pc(i,0)-old));

        old=Pc(i,N-2);
        Pc(i,N-2) = Pc(i,N-2)*(1-as)+
        as*1.0/(  1.0/AUiJ(i,N-2)/*+1.0/AUiJ(i,N-2+1)*//*+1.0/AVIj(i,N-2)+1.0/AVIj(i+1,N-2)  )*(
         U(i,N-2)   +1.0/AUiJ(i,N-2)*   Pc(i,N-2-1)
        //-U(i,N-2+1) +1.0/AUiJ(i,N-2+1)* Pc(i,N-2+1) 
        +V(i,N-2)   +1.0/AVIj(i,N-2)*   Pc(i-1,N-2)
        -V(i+1,N-2) +1.0/AVIj(i+1,N-2)* Pc(i+1,N-2) 
        );
        residual=max(residual,abs(Pc(i,N-2)-old));
}
for (int i=1;i<N-2;i++)
{
    for (int j=1;j<N-2;j++)
    {   
        old=Pc(i,j);
        Pc(i,j) = Pc(i,j)*(1-as)+
        as*1.0/(  1.0/AUiJ(i,j)+1.0/AUiJ(i,j+1)+1.0/AVIj(i,j)+1.0/AVIj(i+1,j)  )*(
         U(i,j)   +1.0/AUiJ(i,j)*   Pc(i,j-1)
        -U(i,j+1) +1.0/AUiJ(i,j+1)* Pc(i,j+1) 
        +V(i,j)   +1.0/AVIj(i,j)*   Pc(i-1,j)
        -V(i+1,j) +1.0/AVIj(i+1,j)* Pc(i+1,j) 
        );
        residual=max(residual,abs(Pc(i,j)-old));
        //0=
        // U(i,j)   +1.0/AUiJ(i,j)*  ( P(i,j-1)-P(i,j) )
        //-U(i,j+1) -1.0/AUiJ(i,j+1)*( P(i,j)-P(i,j+1) )
        //+V(i,j)   +1.0/AVIj(i,j)*  ( P(i-1,j)-P(i,j) )
        //-V(i+1,j) -1.0/AVIj(i+1,j)*( P(i,j)-P(i+1,j) );
    }
}
if (residual < tol/Re/dL)
{
	cout<<"No. Pc loop:  "<<n_Pc_loop<<endl;
	break;
}
}

P=P+ap*Pc;
for (int i=0;i<N-1;i++)
{
    for (int j=0;j<N;j++)
    {
        U(i,j)=U(i,j)+(     (j==0 || j==N-1)?0:( au*1.0/AUiJ(i,j)*(Pc(i,j-1)-Pc(i,j)) )    );
    }
}
for (int i = 0; i < N ; i++)
{
    for (int j = 0; j < N-1; j++)
    {
        V(i, j) = V(i, j) +(    (i==0 || i==N-1)?0:( av * 1.0 / AVIj(i, j) * (Pc(i-1, j) - Pc(i, j)) )    );
    }
}
*/
