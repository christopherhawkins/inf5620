%% function to calculate elements of b
function [b] = calc_b(f, psi, omega)
syms x;
syms X;

h=omega(2)-omega(1);
det_J=h/2;
f_local=subs(f, (omega(2)+omega(1))/2 + X*h/2);

b=int(f_local*psi*det_J, -1, 1);