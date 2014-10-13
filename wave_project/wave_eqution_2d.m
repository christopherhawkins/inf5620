clear;

%%
% simulation variables
nx=100;
ny=100;
dx=0.1;
dy=0.1;
dt=0.05;
end_time=10;
b=0;

max_t=end_time/dt-1;


%%
% Make geometry (array q(x, y)) and unitial velocity field for import to
% program
make_array('geometry.txt', nx, ny, 'cos 2d', round(nx/2), round(ny/2), dx, dy);
make_array('u_initial.txt', nx, ny, 'gaussian 1d', round(nx/2), round(ny/4), dx, dy);


%%
% run main program with variables and output q(x, y), u(x, y). Import data
% as data
[u_raw, u_d_raw]=command(nx, ny, dx, dy, dt, end_time, b, 0);


%% 
% visualisation  set up figure
scrsz = get(0,'ScreenSize');
figure('Position', [500 500 700 500])

h = uicontrol('Style','slider','Position', [20 5 500 20],'Max',max_t,'Min',0,'SliderStep',[1/max_t,10/max_t],'Value',0);
while true
    
     t = round (get(h,'value'));
     
%%
% sort data into 2D array
     for i=1:nx-1
         for j=1:ny-1
             u(j,i)=u_raw(j + i*ny + t*nx*ny); 
             x(j, i)=i;
             y(j, i)=j;
         end
     end
     
%%
% plots go here
     imagesc(u);
     T(t+1)=u(1, 1);
     %plot(x(50, :), u(50, :));
     
     figure(1);
end