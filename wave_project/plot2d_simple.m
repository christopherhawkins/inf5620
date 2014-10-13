%clear


max_t=299;
nx=100;
ny=100;

scrsz = get(0,'ScreenSize');
figure('Position', [500 500 700 500])

h = uicontrol('Style','slider','Position', [20 5 500 20],'Max',max_t,'Min',0,'SliderStep',[1/max_t,10/max_t],'Value',0);
while true
    
     t = round (get(h,'value'));
     for i=1:nx-1
         for j=1:ny-1
             u(j,i)=data(j + i*ny + t*nx*ny); 

             x(j, i)=i;
             y(j, i)=j;
         end
     end
     
     imagesc(u);
     %hold on;
     %plot(u1(:, 50), 'r')
     
     %axis([0 100 -0.1 0.1])
     
     figure(1);
end