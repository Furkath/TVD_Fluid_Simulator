std::clock_t c_start = std::clock();
//std::clock_t c_end = std::clock();

//long_double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
//std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

while (it_num <= max_it && max_res > tol)
{
    it_num++;
    #include "solveMomentumEq_new.H"
    u_old=u;
    v_old=v;
    //cout<<u<<endl;
    //cout<<v<<endl;
    //cout<<d_u<<endl;
    //cout<<d_v<<endl;
    //cout<<"here 2"<<endl;
    #include "solvePressureCorrectEq.H"
    //cout<<"here 3"<<endl;
    //cout<<p<<endl;
    //cout<<p_prime<<endl;
    #include "UV_update.H"
    //cout<<"here 4"<<endl;
    p_star=p;
    v_res=(v-v_old).cwiseAbs().maxCoeff();
    u_res=(u-u_old).cwiseAbs().maxCoeff();
    max_res=max(u_res,v_res);

    cout<<"SIMPLE loop num:  "<<it_num<<" ;   residual:  "<<max_res<<endl;
    //if(it_num>16500 && max_res <0.00005)
    //{
    //	cout<<"将就------"<<endl;
    // 	break;
    //}
    if (max_res > 1.648)
    {
        cout<<"---寄！---"<<endl;
        break;
    }
}


std::clock_t c_end = std::clock();
long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
std::cout << "CPU time used: " << time_elapsed_ms << " ms"<<endl;












/*
//NEW_SIMPLE_LOOP=true;
n_SIMPLE_loop=0;
while (true)
{
    NEW_SIMPLE_LOOP=true;
    n_SIMPLE_loop++;
    //if(is_steady) Info<<"No. SIMPLE loop:  "<<n_SIMPLE_loop<<endl;
    //store the initial or previous loop's result to tell wether converged
    U_ = U;
    V_ = V;
    P_ = P;
    #include "solveMomentumEq.H"

    //cout<<"Between Utop:  "<<U(N-2,20)<<endl;
    //cout<<"After Utop:  "<<U(N-2,0)<<endl;

    #include "solvePressureCorrectEq.H"

    //cout<<"After Utop:  "<<U(N-2,20)<<endl;
    //cout<<"After Utop:  "<<U(N-2,0)<<endl;
    cout<<U<<endl;
    cout<<V<<endl;
    cout<<P<<endl;
    double residual=std::max( (U_ - U).cwiseAbs().maxCoeff(), (V_ - V).cwiseAbs().maxCoeff() );
    residual =std::max( residual , (P_ - P).cwiseAbs().maxCoeff() / Re / dL );
    if(is_steady) cout<<"No. SIMPLE loop:  "<<n_SIMPLE_loop<<endl;
    if(residual<tol)
    {
	    if(!is_steady) cout<<"No. SIMPLE loop:  "<<n_SIMPLE_loop<<endl;
    	    break;
    }
}
*/
