clear;

%%
% simulation variables
nx=100;
ny=100;
dx=0.1;
dy=0.1;
dt=0.01;
end_time=dt*2;
b=0;
A=1;

max_t=end_time/dt-1;


%%
% Make geometry (array q(x, y)) and unitial velocity field for import to
% program
make_array('geometry.txt', nx, ny, 'const 2d', round(nx/2), round(ny/2), dx, dy);
make_array('u_initial.txt', nx, ny, 'manufactured', round(nx/2), round(ny/2), dx, dy);


%%
% run main program with variables and output q(x, y), u(x, y). Import data
% as data
[u_raw, u_d_raw]=command(nx, ny, dx, dy, dt, end_time, b, 2);


c2=1/100;
c1=(-2*c2)/(3*nx*dx);

%%
% sort data into 2D array
for i=1:nx-1
     for j=1:ny-1
         u(j,i)=u_d_raw(j + i*ny); 
         x(j, i)=i;
         y(j, i)=j;
     end
end

%%
% actual derivetive of cubic waveform
for i=1:nx-1
    for j=1:ny-1
        X(i, j)=(1-0.01)*(6*c1*(i*dx)+2*c2)*(c1*((j*dx)^3)+c2*((j*dx)^2));
        Y(i, j)=(1-0.01)*(6*c1*(j*dx)+2*c2)*(c1*((i*dx)^3)+c2*((i*dx)^2));
        a(i, j)=A*(X(i, j)+Y(i, j));
    end
end

%%
% plot derivetives
plot(u(:, 20));
hold on;
plot(a(:, 21), 'r');