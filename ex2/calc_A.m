%% function to calculate elements of A 
function [A] = calc_A(psi_i, psi_j, omega)

h=omega(2)-omega(1);
det_J=h/2;
A=int(psi_i*psi_j*det_J, -1, 1);