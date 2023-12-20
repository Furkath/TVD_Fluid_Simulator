if(TVD && ii==jj)
{
    int N=ii-1;
    // Construct surface fluxes
    MatrixXd FUu(N - 1, N + 1);
    MatrixXd FUv(N, N);
    MatrixXd FVu(N, N);
    MatrixXd FVv(N + 1, N - 1);
    // Construct r in the TVD flux limiter function for U
    MatrixXd rUIsJ_plus(N, N);
    MatrixXd rUIsJ_minus(N, N);
    MatrixXd rUijs_plus(N - 1, N + 1);
    MatrixXd rUijs_minus(N - 1, N + 1);
    // Construct the souurce term for U
    MatrixXd alphaUu(N - 1, N + 1);
    MatrixXd alphaUv(N, N);
    MatrixXd SUiJ(N - 1, N);
    // Construct r in the TVD flux limiter function for V
    MatrixXd rVisj_plus(N + 1, N - 1);
    MatrixXd rVisj_minus(N + 1, N - 1);
    MatrixXd rVIJs_plus(N, N);
    MatrixXd rVIJs_minus(N, N);
    // Construct the source term for V
    MatrixXd alphaVv(N + 1, N - 1);
    MatrixXd alphaVu(N, N);
    MatrixXd SVIj(N, N - 1);

    int u_row=u.rows()-2;
    int u_col=u.cols()-2;
    int v_row=v.rows()-2;
    int v_col=v.cols()-2;
    MatrixXd Utmp = u.block(1,1,u_row,u_col);
    MatrixXd U=Utmp.transpose();
    MatrixXd Vtmp = v.block(1,1,v_row,v_col);
    MatrixXd V=Vtmp.transpose();
    //Re=rho*1.0*1.0/mu;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            if (j == 0 || j == N)
                { FUu(i, j) = 0; }// B.C.
            else
                { FUu(i, j) = 0.5 * (U(i, j - 1) + U(i, j)); }
        }
    }
    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            if (i == 0 || i == N)
                { FVv(i, j) = 0; }// B.C.
            else
                { FVv(i, j) = 0.5 * (V(i - 1, j) + V(i, j)); }
        }
    }
    for (int i = 0; i < N; i++)
    {
        FUv(i, 0) = 0.5 * (0 + V(i, 0)); // B.C.
        for (int j = 1; j < N - 1; j++)
        {
            FUv(i, j) = 0.5 * (V(i, j - 1) + V(i, j));
        }
        FUv(i, N - 1) = 0.5 * (V(i, N - 2) + 0); // B.C.
    }
    for (int j = 0; j < N; j++)
    {
        FVu(0, j) = 0.5 * (0 + U(0, j)); // B.C.
        for (int i = 1; i < N - 1; i++)
        {
            FVu(i, j) = 0.5 * (U(i - 1, j) + U(i, j));
        }
        FVu(N - 1, j) = 0.5 * (U(N - 2, j) + 1.0); // 1.0 is the upper lid velocity as normalized
    }

    // update alpha_U and alpha_V
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            alphaUv(i, j) = ((FUv(i, j) >= 0) ? 1 : 0);
            alphaVu(i, j) = ((FVu(i, j) >= 0) ? 1 : 0);
        }
    }
    for (int i = 0; i < N - 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            alphaUu(i, j) = ((FUu(i, j) >= 0) ? 1 : 0);
        }
    }
    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 0; j < N - 1; j++)
        {
            alphaVv(i, j) = ((FVv(i, j) >= 0) ? 1 : 0);
        }
    }

    // update the TVD r for U
    for (int j = 0; j < N; j++)
    {
        rUIsJ_plus(0, j) = 0;                                           // Consider the B.C.
        rUIsJ_minus(0, j) = (U(0 + 1, j) - U(0, j)) / (U(0, j) - 0);    // B.C.
        rUIsJ_plus(1, j) = (0 - U(1 - 1, j)) / (U(1 - 1, j) - U(1, j)); // B.C.
        rUIsJ_minus(1, j) = (U(1 + 1, j) - U(1, j)) / (U(1, j) - U(1 - 1, j));
        for (int i = 2; i < N - 2; i++)
        {
            rUIsJ_plus(i, j) = (U(i - 2, j) - U(i - 1, j)) / (U(i - 1, j) - U(i, j));
            rUIsJ_minus(i, j) = (U(i + 1, j) - U(i, j)) / (U(i, j) - U(i - 1, j));
        }
        rUIsJ_plus(N - 2, j) = (U(N - 2 - 2, j) - U(N - 2 - 1, j)) / (U(N - 2 - 1, j) - U(N - 2, j));
        rUIsJ_minus(N - 2, j) = (1.0 - U(N - 2, j)) / (U(N - 2, j) - U(N - 2 - 1, j));        // Consider the B.C.
        rUIsJ_plus(N - 1, j) = (U(N - 1 - 2, j) - U(N - 1 - 1, j)) / (U(N - 1 - 1, j) - 1.0); // B.C.
        rUIsJ_minus(N - 1, j) = 0;                                                            // B.C.
    }
    for (int i = 0; i < N - 1; i++)
    {
        rUijs_plus(i, 0) = 0;                                                         // Consider the B.C.
        rUijs_minus(i, 0) = 0;                                                        // B.C.
        rUijs_plus(i, 1) = (2 * 0 - U(i, 0) - U(i, 1 - 1)) / (U(i, 1 - 1) - U(i, 1)); // B.C. type1
        rUijs_minus(i, 1) = (U(i, 1 + 1) - U(i, 1)) / (U(i, 1) - U(i, 1 - 1));
        for (int j = 2; j < N - 1; j++)
        {
            rUijs_plus(i, j) = (U(i, j - 2) - U(i, j - 1)) / (U(i, j - 1) - U(i, j));
            rUijs_minus(i, j) = (U(i, j + 1) - U(i, j)) / (U(i, j) - U(i, j - 1));
        }
        rUijs_plus(i, N - 1) = (U(i, N - 1 - 2) - U(i, N - 1 - 1)) / (U(i, N - 1 - 1) - U(i, N - 1));
        rUijs_minus(i, N - 1) = (2 * 0 - U(i, N - 1) - U(i, N - 1)) / (U(i, N - 1) - U(i, N - 1 - 1)); // Consider the B.C. type1
        rUijs_plus(i, N) = 0;                                                                          // B.C.
        rUijs_minus(i, N) = 0;                                                                         // B.C.
    }
    // treat B.C. directly

    // update the TVD r for V
    for (int j=0;j<N-1;j++)
    {
        rVisj_plus(0, j) = 0;                                                         // Consider the B.C.
        rVisj_minus(0, j) = 0;                                                        // B.C.
        rVisj_plus(1, j) = (2 * 0 - V(0, j) - V(1 - 1, j)) / (V(1 - 1, j) - V(1, j)); // B.C. type1
        rVisj_minus(1, j) = (V(1 + 1, j) - V(1, j)) / (V(1, j) - V(1 - 1, j));
        for (int i = 2; i < N - 1; i++)
        {
            rVisj_plus(i, j) = (V(i - 2, j) - V(i - 1, j)) / (V(i - 1, j) - V(i, j));
            rVisj_minus(i, j) = (V(i + 1, j) - V(i, j)) / (V(i, j) - V(i - 1, j));
        }
        rVisj_plus(N - 1, j) = (V(N - 1 - 2, j) - V(N - 1 - 1, j)) / (V(N - 1 - 1, j) - V(N - 1, j));
        rVisj_minus(N - 1, j) = (2 * 0 - V(N - 1, j) - V(N - 1, j)) / (V(N - 1, j) - V(N - 1 - 1, j)); // Consider the B.C. type1
        rVisj_plus(N, j) = 0;                                                                          // B.C.
        rVisj_minus(N, j) = 0;                                                                         // B.C.
    }
    for (int i = 0; i < N; i++)
    {
        rVIJs_plus(i, 0) = 0;                                           // B.C.
        rVIJs_minus(i, 0) = (V(i, 0 + 1) - V(i, 0)) / (V(i, 0) - 0);    // B.C.
        rVIJs_plus(i, 1) = (0 - V(i, 1 - 1)) / (V(i, 1 - 1) - V(i, 1)); // B.C.
        rVIJs_minus(i, 1) = (V(i, 1 + 1) - V(i, 1)) / (V(i, 1) - V(i, 1 - 1));
        for (int j = 2; j < N - 2; j++)
        {
            rVIJs_plus(i, j) = (V(i, j - 2) - V(i, j - 1)) / (V(i, j - 1) - V(i, j));
            rVIJs_minus(i, j) = (V(i, j + 1) - V(i, j)) / (V(i, j) - V(i, j - 1));
        }
        rVIJs_plus(i, N - 2) = (V(i, N - 2 - 2) - V(i, N - 2 - 1)) / (V(i, N - 2 - 1) - V(i, N - 2));
        rVIJs_minus(i, N - 2) = (0 - V(i, N - 2)) / (V(i, N - 2) - V(i, N - 2 - 1));        // B.C.
        rVIJs_plus(i, N - 1) = (V(i, N - 1 - 2) - V(i, N - 1 - 1)) / (V(i, N - 1 - 1) - 0); // B.C.
        rVIJs_minus(i, N - 1) = 0;                                                          // B.C.
    }
    // treat B.C. directly

    // calculate Sij
    for (int i=0;i<N-1;i++)
    {
        for (int j=0;j<N;j++)
        {
            SUiJ(i,j)=0.5*(                                                                                                                                   
                      FUv(i+1,j)*(   (1-alphaUv(i+1,j)) * psi(rUIsJ_minus(i+1,j))  -  alphaUv(i+1,j)*psi(rUIsJ_plus(i+1,j))   )*(  (i==N-2)  ?   (  1.0-U( i,j )  ):(  U( i+1,j )-U( i,j )  )   )    
                     +FUv(i,j)*(    alphaUv(i,j) * psi(rUIsJ_plus(i,j)) - (1-alphaUv(i,j)) * psi(rUIsJ_minus(i,j))   )*(           (i==0)    ?   (  U( i,j )-0.0  ):(  U( i,j )-U( i-1,j )  )   )                                                                                                                       
                     +FUu(i,j+1)*(   (1-alphaUu(i,j+1)) * psi(rUijs_minus(i,j+1))  -  alphaUu(i,j+1)*psi(rUijs_plus(i,j+1))   )*(  (j==N-1)  ?   (  2*0.0-2*U( i,j )  ):(  U( i,j+1 )-U( i,j )  )   )   
                     +FUu(i,j)*(    alphaUu(i,j) * psi(rUijs_plus(i,j)) - (1-alphaUu(i,j)) * psi(rUijs_minus(i,j))   )*(           (j==0 )   ?   (  2*U( i,j )-2*0.0  ):(  U( i,j )-U( i,j-1 )  )   )                                                                                                                                             
            );
        }
    }
    for (int i = 0; i < N ; i++)
    {
        for (int j = 0; j < N-1; j++)
        {
            SVIj(i, j) =0.5*(                                                                                                                                   
                      FVv(i+1,j)*(   (1-alphaVv(i+1,j)) * psi(rVisj_minus(i+1,j))  -  alphaVv(i+1,j)*psi(rVisj_plus(i+1,j))   )*(  (i==N-1)  ?   (  2*0.0-2*V( i,j )  ):(  V( i+1,j )-V( i,j )  )   )     
                     +FVv(i,j)*(    alphaVv(i,j) * psi(rVisj_plus(i,j)) - (1-alphaVv(i,j)) * psi(rVisj_minus(i,j))   )*(           (i==0)    ?   (  2*V( i,j )-2*0.0  ):(  V( i,j )-V( i-1,j )  )   )                                                                                                                              
                     +FVu(i,j+1)*(   (1-alphaVu(i,j+1)) * psi(rVIJs_minus(i,j+1))  -  alphaVu(i,j+1)*psi(rVIJs_plus(i,j+1))   )*(  (j==N-2)  ?   (  0.0-V( i,j )  ):(  V( i,j+1 )-V( i,j )  )   )      
                     +FVu(i,j)*(    alphaVu(i,j) * psi(rVIJs_plus(i,j)) - (1-alphaVu(i,j)) * psi(rVIJs_minus(i,j))   )*(           (j==0)    ?   (  V( i,j )-0.0  ):(  V( i,j )-V( i,j-1 )  )   )                                                                                                                                                  
            );
        }
    }

    add_su = SUiJ.transpose();
    add_sv = SVIj.transpose();

}
