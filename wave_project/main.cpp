#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>

using namespace std;

//--------------------------------------------------------------------------
// global variables
double dx, dy, dt;
int nx, ny;

//--------------------------------------------------------------------------
// function to acquire parameters
double get_param_double(int argc, char *argv[], const char *b){
	int i;
	double var;
	bool equal;

	for (i=0;i<argc;i++){
		string c1 = argv[i], c2 = b;
		equal = (c1 == c2);
		if(equal){
			var=atof(argv[i+1]);
		}
	}

	return var;
}

//--------------------------------------------------------------------------
// function to import file to 2D array
void import_file_2D(ifstream *file, double **a){
	int i, j;

	if(file->is_open()){
		for(i=0; i<nx; i++)
			for(j=0; j<ny; j++){
				file->operator >>(a[i][j]);
			}
	}
}

//---------------------------------------------------------------------------
// spacial derivative u''' in x
double d3_dx3(double **u, int i, int j){
	return ((-1/2)*u[i-2][j] + u[i-1][j] - u[i+1][j] + (1/2)*u[i+2][j])/(dx*dx*dx);
}

//---------------------------------------------------------------------------
// spacial derivative u''' in y
double d3_dy3(double **u, int i, int j){
	return ((-1/2)*u[i][j-2] + u[i][j-1] - u[i][j+1] + (1/2)*u[i][j+2])/(dy*dy*dy);
}

//---------------------------------------------------------------------------
// spacial derivative u'' in x
double d2_dx2(double **u, int i, int j){
	return (u[i+1][j] - 2*u[i][j] + u[i-1][j])/(dx*dx);
}

//--------------------------------------------------------------------------
// spacial derivative u'' in y
double d2_dy2(double **u, int i, int j){
	return (u[i][j+1] - 2*u[i][j] + u[i][j-1])/(dy*dy);
}

//--------------------------------------------------------------------------
// spacial derivative u' in x
double d_dx(double **u, int i, int j){
	return (u[i+1][j]-u[i][j])/(dx);
}

//--------------------------------------------------------------------------
// spacial derivative u' in y
double d_dy(double **u, int i, int j){
	return (u[i][j+1]-u[i][j])/(dy);
}

//--------------------------------------------------------------------------
// isotropic spacial derivative u' in x
double id_dx(double **u, int i, int j){
	return (u[i+1][j+1] + 2*u[i+1][j] - u[i-1][j-1] - u[i-1][j+1] - 2*u[i-1][j] + u[i+1][j-1])/(8*dx);
}

//--------------------------------------------------------------------------
// isotropic spacial derivative u' in y
double id_dy(double **u, int i, int j){
	return (u[i+1][j+1] + 2*u[i][j+1] - u[i-1][j-1] + u[i-1][j+1] - 2*u[i][j-1] - u[i+1][j-1])/(8*dy);
}

//--------------------------------------------------------------------------
// determine spacial derivative of equation with finite difference. chain rule expanded
void grad_dot_q_grad_c(double **u_d, double **u, double **q){
	int i, j;

	for(i=1; i<nx-1; i++)
		for(j=1; j<ny-1; j++)
			u_d[i][j]=q[i][j]*(d2_dx2(u, i, j) + d2_dy2(u, i, j)) + id_dx(q, i, j)*id_dx(u, i, j) + id_dy(q, i, j)*id_dy(u, i, j);
}

//--------------------------------------------------------------------------
// determine spacial derivative of equation with finite difference.
void grad_dot_q_grad(double **u_d, double **u, double **q){
	int i, j;
	double q_ph_x, q_mh_x, q_ph_y, q_mh_y;

	for(i=1; i<nx-1; i++)
		for(j=1; j<ny-1; j++){
			q_ph_x=(q[i+1][j] + q[i][j])/2;			x2=x*x;
			x3=x*x*x;
			y=(j+1)*dy;
			y2=y*y;
			y3=y*y*y;
			q_mh_x=(q[i-1][j] + q[i][j])/2;

			q_ph_y=(q[i][j+1] + q[i][j])/2;
			q_mh_y=(q[i][j-1] + q[i][j])/2;

			u_d[i][j]=(q_ph_x*(u[i+1][j] - u[i][j]) - q_mh_x*(u[i][j] - u[i-1][j]))/(dx*dx) + (q_ph_y*(u[i][j+1] - u[i][j]) - q_mh_y*(u[i][j] - u[i][j-1]))/(dy*dy);
		}
}

//--------------------------------------------------------------------------
// time integration via standard finite difference
void integrate_standard(double **u_n, double **u, double **u_d, double **u_t, double b){
	int i, j;
	double dt2=dt*dt;
	double dt_b=dt*b;

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
			u_n[i][j]=((dt2)*u_d[i][j] + (2-dt_b)*u[i][j] - u_t[i][j])/(1+dt_b);
}

//--------------------------------------------------------------------------
// boundary conditions
void boundary_zero_gradient(double **u_d){
	int i, j;

	for(i=0; i<nx; i++){
		u_d[i][0]=u_d[i][1];
		u_d[i][ny-1]=u_d[i][ny-2];
	}

	for(j=0; j<ny; j++){
		u_d[0][j]=u_d[1][j];
		u_d[nx-1][j]=u_d[nx-2][j];
	}
}

//--------------------------------------------------------------------------
// update time vectors
void update(double **u_n, double **u, double **u_t){
	int i, j;

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++){
			u_t[i][j]=u[i][j];
			u[i][j]=u_n[i][j];
		}
}

//--------------------------------------------------------------------------
// manufactured solution
void manufactured_solution(double **u, int t){
	int i, j;
	double x, L;


	A=1;
	c2=1./100;
	c1=(-2.*c2)/(3*nx*dx);

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++){
			x=(i+1)*dx;
			L=nx*dx;
			u[i][j]+=(2-x*(L-x))*sin(t);
		}
}

//--------------------------------------------------------------------------
// cubic test source term
void cubic_test(double **u, int t){
	int i, j;
	double c2, c1;
	double x, x2, x3, y, y2, y3;
	double A;

	A=1;
	c2=1./100;
	c1=(-2.*c2)/(3*nx*dx);

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++){
			x=(i+1)*dx;
			x2=x*x;
			x3=x*x*x;
			y=(j+1)*dy;
			y2=y*y;
			y3=y*y*y;

			u[i][j]+=A*(dt*t-1)*((c1*y3+c2*y2)*(c1*6*x + c2*2) + (c1*x3+c2*x2)*(c1*6*y + c2*2));
		}
}

//--------------------------------------------------------------------------
// truncation error
void truncation_error(double **t_e, double **g_x, double **g_y, double **u, double **q){
	int i, j;

	for(i=2; i<nx-2; i++)
		for(j=2; j<ny-2; j++){
			g_x[i][j]=(1./24)*q[i][j]*d3_dx3(u, i, j)*(dx*dx) + (1./8)*d_dx(u, i, j)*d2_dx2(q, i, j)*(dx*dx);
			g_y[i][j]=(1./24)*q[i][j]*d3_dx3(u, i, j)*(dx*dx) + (1./8)*d_dx(u, i, j)*d2_dx2(q, i, j)*(dx*dx);
		}

	for(i=2; i<nx-2; i++)
		for(j=2; j<ny-2; j++){
			t_e[i][j]=d_dx(g_x, i, j)*dx*dx+d_dy(g_y, i, j)*dy*dy;
		}
}

//--------------------------------------------------------------------------
// main program
int main(int argc, char *argv[]){

//--------------------------------------------------------------------------
// initiate variables
	double **u, **u_t, **u_d, **u_n;
	double **q;
	double **g_x, **g_y, **t_e;
	double b;
	double end_time;
	double test;

	int i, j;
	int t;

	bool running;

	b=get_param_double(argc, argv, "-b");
	nx=get_param_double(argc, argv, "-nx");
	ny=get_param_double(argc, argv, "-ny");
	dx=get_param_double(argc, argv, "-dx");
	dy=get_param_double(argc, argv, "-dy");
	dt=get_param_double(argc, argv, "-dt");
	end_time=get_param_double(argc, argv, "-end_time");
	test=get_param_double(argc, argv, "-test");

	u=new double *[nx];
	u_t=new double *[nx];
	u_d=new double *[nx];
	u_n=new double *[nx];

	g_x=new double *[nx];
	g_y=new double *[nx];
	t_e=new double *[nx];

	q=new double *[nx];

	for(i=0; i<nx; i++){
		u[i]=new double [ny];
		u_t[i]=new double [ny];
		u_d[i]=new double [ny];
		u_n[i]=new double [ny];

		g_x[i]=new double [ny];
		g_y[i]=new double [ny];
		t_e[i]=new double [ny];

		q[i]=new double [ny];

	}

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++){
			u[i][j]=0;
			u_t[i][j]=0;
			u_d[i][j]=0;
			u_n[i][j]=0;
			q[i][j]=0;
			g_x[i][j]=0;
			g_y[i][j]=0;
			t_e[i][j]=0;
		}

	t=0;

//--------------------------------------------------------------------------
// setup the output files
	ofstream u_out("u.txt");
	ofstream u_d_out("u_d.txt");

//--------------------------------------------------------------------------
// import geometry and initial conditions from file
	ifstream geometry("geometry.txt");
	ifstream u_initial("u_initial.txt");

	import_file_2D(&geometry, q);
	import_file_2D(&u_initial, u);

	geometry.close();
	u_initial.close();

	for(i=0; i<nx; i++)
		for(j=0; j<ny; j++){
			u_t[i][j]=u[i][j];
		}

//--------------------------------------------------------------------------
// set main loop running
	running=true;
	while(running){

//--------------------------------------------------------------------------
// calculate the spacial derivative
		grad_dot_q_grad_c(u_d, u, q);

//--------------------------------------------------------------------------
// boundary conditions
		boundary_zero_gradient(u_d);

//--------------------------------------------------------------------------
// tests and errors
		if(test==1){
			cubic_test(u_d, t);
			truncation_error(t_e, g_x, g_y, u, q);
		}
		
		if(test==2){
		        manufactured_solution(u_d, t);
			truncation_error(t_e, g_x, g_y, u, q);
		}
		

//--------------------------------------------------------------------------
// integrate in time
		integrate_standard(u_n, u, u_d, u_t, b);

//--------------------------------------------------------------------------
// update old time arrays
		update(u_n, u, u_t);

//--------------------------------------------------------------------------
// write data to file
		for(i=0; i<nx; i++)
			for(j=0; j<ny; j++){
				u_out << u[i][j] << "\n";
				u_d_out << u_d[i][j] << "\n";
			}

//--------------------------------------------------------------------------
// iterate time in time steps
		t+=1;

//--------------------------------------------------------------------------
// check for end time and stop running if reached
		if(t*dt>=end_time){
			running=false;
		}
	}

//--------------------------------------------------------------------------
// close files
	u_out.close();
	u_d_out.close();

//--------------------------------------------------------------------------
// terminate program
	return 0;
}
