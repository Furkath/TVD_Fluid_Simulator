#include "readinput.H"

int it_num =0;

double mu=rho*v_top*L/Re;
const double dx=L/(ii-1);
const double dy=L/(jj-1);


MatrixXd p = MatrixXd::Zero(ii, jj);
MatrixXd p_star = MatrixXd::Zero(ii, jj);
MatrixXd p_prime = MatrixXd::Zero(ii, jj);
MatrixXd bp = MatrixXd::Zero(ii, jj);
//MatrixXd div = MatrixXd::Zero(ii, jj);

MatrixXd v_star = MatrixXd::Zero(ii, jj + 1);
MatrixXd v_old = MatrixXd::Zero(ii, jj + 1);
double v_res;
MatrixXd v = MatrixXd::Zero(ii, jj + 1);
MatrixXd d_v = MatrixXd::Zero(ii, jj + 1);

MatrixXd u_star = MatrixXd::Zero(ii + 1, jj);
MatrixXd u_old = MatrixXd::Zero(ii + 1, jj);
double u_res;
MatrixXd u = MatrixXd::Zero(ii + 1, jj);
MatrixXd d_u = MatrixXd::Zero(ii + 1, jj);

MatrixXd u_0 = MatrixXd::Zero(ii + 1, jj);
MatrixXd v_0 = MatrixXd::Zero(ii, jj + 1);
//const double endt;
double dt;

for (int i=0;i<ii+1;i++)
{
    u_star(i,jj-1)=v_top;
    u(i,jj-1)=v_top;
}



/*
// Set grid size
const int N = 5;
// Calculate the cell length
double dL = 1.0 / N;
// Set Reynold number
double Re=1;
// calculate the actual (non-normalized) velocity of the lid
double V_actual = Re / 1000.0; //[m/s]
// set the time step
double dt;
// set the total time
double time;
//Define an identifier to tell wether this is the stage for updating fluxes at every new SIMPLE loop 
bool NEW_SIMPLE_LOOP;
//Define an identifier to update U, V seperately
bool UorV;//U1V0
//Set the relaxation parameter for Gauss-Sieldel iterration
double as=0.9;
// Set the relaxation parameter of U for transient SIMPLE loop 
double au=0.8;
// Set the relaxation parameter of V for transient SIMPLE loop
double av = 0.8;
// Set the relaxation parameter of P for transient SIMPLE loop
double ap = 0.8;
// A steady state problem or transient
// bool is_transient;
// Define loop nuimber
int n_SIMPLE_loop;
int n_U_loop;
int n_V_loop;
int n_Pc_loop;
//set residual
double tol=1e-6;
// Construct historical storage, or just write into files, for results after every wanted time interval
// std::vector<Eigen::Matrix<double, N - 1, N - 1>> P_t;
// std::vector<Eigen::Matrix<double, N - 1, N>> U_t;
// std::vector<Eigen::Matrix<double, N, N - 1>> V_t;
// Construct unknown fields
MatrixXd P(N - 1, N - 1);
MatrixXd Pc(N-1,N-1);
MatrixXd U(N - 1, N);
MatrixXd V(N, N - 1);
// Store previous time step's results
MatrixXd U_0(N - 1, N);
MatrixXd V_0(N, N - 1);
// Store previous SIMPLE ieration's results
MatrixXd U_(N - 1, N);
MatrixXd V_(N, N - 1);
MatrixXd P_(N-1,N-1);
// Construct surface fluxes
MatrixXd FUu(N - 1, N + 1);
MatrixXd FUv(N, N);
MatrixXd FVu(N, N);
MatrixXd FVv(N + 1, N - 1);
// Construct central and neighbour coefficients for U
MatrixXd AUi_1J(N - 1, N);
MatrixXd AUi1J(N - 1, N);
MatrixXd AUiJ_1(N - 1, N);
MatrixXd AUiJ1(N - 1, N);
MatrixXd AUiJ(N - 1, N);
// Construct r in the TVD flux limiter function for U
MatrixXd rUIsJ_plus(N, N);
MatrixXd rUIsJ_minus(N, N);
MatrixXd rUijs_plus(N - 1, N + 1);
MatrixXd rUijs_minus(N - 1, N + 1);
// Construct the souurce term for U
MatrixXd alphaUu(N - 1, N + 1);
MatrixXd alphaUv(N, N);
MatrixXd SUiJ(N - 1, N);
// Construct central and neighbour coefficients for V
MatrixXd AVI_1j(N, N - 1);
MatrixXd AVI1j(N, N - 1);
MatrixXd AVIj_1(N, N - 1);
MatrixXd AVIj1(N, N - 1);
MatrixXd AVIj(N, N - 1);
// Construct r in the TVD flux limiter function for V
MatrixXd rVisj_plus(N + 1, N - 1);
MatrixXd rVisj_minus(N + 1, N - 1);
MatrixXd rVIJs_plus(N, N);
MatrixXd rVIJs_minus(N, N);
// Construct the source term for V
MatrixXd alphaVv(N + 1, N - 1);
MatrixXd alphaVu(N, N);
MatrixXd SVIj(N, N - 1);
*/
