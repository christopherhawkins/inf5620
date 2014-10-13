q=load('geometry.txt');
load('mycmap');

surf(u, 'edgecolor', 'none');
camlight left;
view(-58,36);
%axis([0 100 0 100 0 0.23])
alpha(.6);
set(gcf,'Colormap', mycmap)

hold on;

surf(q-0.4, 'edgecolor', 'none');