   ifstream fin("Parameters");
   std::string line;
   std::vector<double> par;
   if(fin)
   {
       while(getline(fin,line))
       {
           //if(line.size()==0)
           //{
           //    break;
           //}
	   std::string ifend;
	   stringstream testend;
   	   testend<<line;
	   testend>>ifend;
	   if(ifend=="END")
	   {
		   break;
	   }
	   else if(ifend!="")
	   {
           double un;
           stringstream sun;
           sun << line;
           sun >> un;
	   	if(un!=0.0)
	   	{
           	par.push_back(un);
	   	}
	   }
       }
   }
   else
   {
       cout<<"no such file"<<endl;
   }
   fin.close();

   //for(int i=0;i<u.size();i++)
   //{
   //     cout<<u[i]<<endl;
   //}
   
//const bool TVD=false;
//const int ii=160;
//const int jj=160;
const int max_it=par[0];	//180000;
double max_res=par[1];	//1000.0;
double tol=par[2];		//1e-05;//1.25575e-05;//2.53e-05;//1.001e-05;//if S 0.8

double alphaP=par[3];	//0.095;//0.095;//0.1;//0.1;//0.01
double alphaS=par[4];	//0.9;//0.8;//0.7;//0.2;//0.02
double alphaU=par[5];	//0.3;//0.95;//0.8;//0.3;//0.03
double alphaV=par[6];	//0.3;//0.9;//0.8;//0.3;//0.03

double rho=par[7];		//1.0;
double v_top=par[8];		//1.0;
double L=par[9];		//1.0;

double Re=par[10];
const int ii=par[11];		//160;
const int jj=par[12];		//160;
const bool TVD=(par[13]>0?true:false);		//false;
const bool SorT=(par[14]>0?true:false);
const double endt=par[15];
double writet=par[16];

if(SorT) writet=-1.0;
if(writet>=endt) writet=endt;

//double Re=100;
        //double mu=0.01;
//double rho=1.0;
//double v_top=1.0;
//double L=1.0;



//double alphaP=0.095;//0.095;//0.1;//0.1;//0.01
//double alphaS=0.9;//0.8;//0.7;//0.2;//0.02
//double alphaU=0.3;//0.95;//0.8;//0.3;//0.03
//double alphaV=0.3;//0.9;//0.8;//0.3;//0.03
        //double alphas=0.625*alphaS;

//double tol=1e-05;//1.25575e-05;//2.53e-05;//1.001e-05;//if S 0.8
 

