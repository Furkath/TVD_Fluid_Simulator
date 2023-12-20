{
v = Eigen::MatrixXd::Zero(ii, jj + 1);
u = Eigen::MatrixXd::Zero(ii + 1, jj);
for (int i=1;i<ii;i++)
{
    for (int j=1;j<jj-1;j++)
    {
        u(i,j) = u_star(i,j) + alphaU*d_u(i,j)*(p_prime(i-1,j)-p_prime(i,j));
    }
}
for  (int i=1;i<ii-1;i++)
{
    for (int j=1;j<jj;j++)
    {
        v(i,j) = v_star(i,j) + alphaV*d_v(i,j)*(p_prime(i,j-1)-p_prime(i,j));
    }
}




v.row(0)=Eigen::VectorXd::Zero(jj+1);
v.row(ii-1)=Eigen::VectorXd::Zero(jj+1);
v.col(0)=-v.col(1);
v.col(jj)=-v.col(jj-1);

u.row(0)=-u.row(1);
u.row(ii)=-u.row(ii-1);
u.col(0)=Eigen::VectorXd::Zero(ii+1);
u.col(jj-1)=v_top*Eigen::VectorXd::Ones(ii+1);

}
