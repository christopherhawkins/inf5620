function[]=make_array(file_name, nx, ny, type, bx, by, dx, dy)

A=0.1;
c=10;
background=0.2;

for i=1:nx
    for j=1:ny
        a(i, j)=0;
    end
end


%%
if(strcmp(type, 'gaussian 1d'))
    for i=1:nx
         for j=1:ny
            a(i, j)=A*exp(-(j-by)^2/(2*c^2));
            a(i, j)=a(i, j)+0*background;
         end
    end
end

%%
if(strcmp(type, 'gaussian 2d'))
    for i=1:nx
         for j=1:ny
            a(i, j)=background;
            a(i, j)=a(i, j)+A*exp(-(j-by)^2/(2*c^2))*exp(-(i-bx)^2/(2*c^2));
         end
    end
end

%%
if(strcmp(type, 'cos 2d'))
    for i=1:nx
         for j=1:ny
             a(i, j)=background;
             if((i-bx)^2 + (j-by)^2 < c^2)
                a(i, j)=a(i, j)+A*(cos(pi*(i-bx)/(2*c)))*(cos(pi*(j-by)/(2*c)));
             end
         end
    end
end

%%
if(strcmp(type, 'sharp 2d'))
    for i=1:nx
         for j=1:ny
             a(i, j)=background;
             if((i-bx)^2 + (j-by)^2 < c^2)
                a(i, j)=a(i, j)+A;
             end
         end
    end
end

%%
if(strcmp(type, 'const 2d'))
    for i=1:nx
         for j=1:ny
             a(i, j)=1;
         end
    end
end

%%
if(strcmp(type, 'test'))
    A=1;
    c2=1/100;
    c1=(-2*c2)/(3*nx*dx);

    for i=1:nx
         for j=1:ny
             X(i, j)=c1*((i*dx)^3)+c2*((i*dx)^2);
             Y(i, j)=c1*((j*dx)^3)+c2*((j*dx)^2);
             a(i, j)=A*X(i, j)*Y(i, j);
         end
    end
end

%%
if(strcmp(type, 'standing'))
    dx=0.1;
    w=pi;
    kx=pi/(nx*dx);
    ky=pi/(ny*dy);

    for i=1:nx
         for j=1:ny
             x=i*dx;
             y=j*dy;
             a(i, j)=A*cos(x*kx)*cos(y*ky);
         end
    end
end

%%
dlmwrite(file_name, a, 'delimiter',' ');