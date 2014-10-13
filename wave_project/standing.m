clear;
%%
% simulation variables
h=0.01:0.01:0.1;
nx=100;
ny=100;
end_time=1;
dx=0.1;
dy=0.1;

A=0.1;
b=0.1;
kx=pi/(nx*dx);
ky=pi/(ny*dy);
w=sqrt(kx^2+ky^2);
w_damped=sqrt(kx^2+ky^2);

%%
% exact solution
for i=1:nx
    for j=1:ny
        x=i*dx;
        y=j*dy;
        u_e(j, i)=A*cos(x*kx)*cos(y*ky)*cos(w*end_time);
        u_e_damped(j, i)=A*exp(-b*end_time)*cos(x*kx)*cos(y*ky)*cos(w_damped*end_time);
    end
end

%%
% create arrays for import
make_array('geometry.txt', nx, ny, 'const 2d', round(nx/2), round(ny/2), dx, dy);
make_array('u_initial.txt', nx, ny, 'standing', round(nx/2), round(ny/2), dx, dy);

%%
% run program for each time spacing and compute error
for t=1:length(h)
    
    [u_raw, u_d_raw]=command(nx, ny, dx, dy, h(t), end_time, b, 0);
    
    err(t)=0;
    
    for i=1:nx-1
        for j=1:ny-1
            u(j,i)=u_raw(j + i*ny + round(end_time/h(t)-1)*nx*ny);
            err(t)=err(t)+sqrt((u(j,i)-u_e(j,i))^2);
            if(b~=0)
                err(t)=err(t)+sqrt((u(j,i)-u_e_damped(j,i))^2);
            end
        end
    end
end

plot(err)