function[u, u_d]=command(nx, ny, dx, dy, dt, end_time, b, test)

cstr1='./inf5620_project2';
cstr2='-b';
cstr3='-nx';
cstr4='-ny';
cstr5='-dx';
cstr6='-dy';
cstr7='-dt';
cstr8='-end_time';
cstr9='-test';

command_string=strcat(cstr1, 32, cstr2, 32, num2str(b), 32, cstr3, 32, num2str(nx), 32, cstr4, 32, num2str(ny), 32, cstr5, 32, num2str(dx), 32, cstr6, 32, num2str(dy), 32, cstr7, 32, num2str(dt), 32, cstr8, 32, num2str(end_time), 32, cstr9, 32, num2str(test));
system(command_string);

u=load('u.txt');
u_d=load('u_d.txt');
