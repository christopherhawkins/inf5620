# -*- coding: utf-8 -*-
import sympy as sp
import numpy as np
#--------------------------------------------------------------------------
def calc_A(psi_i, psi_j, omega):

    h=omega[1] - omega[0];
    det_J=h/2;
    return sp.integrate(psi_i*psi_j*det_J, (X, -1, 1));

#--------------------------------------------------------------------------
def calc_b(f, psi, omega):
    X=sp.Symbol('X');

    h=omega[1] - omega[0];
    det_J=h/2;
    x = (omega[0] + omega[1])/2 + h/2*X;
    f_local = f.subs(’x’, x);
    return sp.integrate(f_local*psi*det_J, (X, -1, 1));

#--------------------------------------------------------------------------
def P(d):
    X=sp.Symbol('X');

    domain=[-1, 1];
    for e in range(0, d):
	Pa[e]=domain[0] + (e-1)*(domain[1]-domain[0])/(d);
    end

    for r in range(0, d):
	psi[r]=0*X+1;
	for s in range(0, d):
	    if r!=s:
		psi[r]=psi[r]*(X-Pa[s])/(Pa[r] - Pa[s]);

    return psi;

#--------------------------------------------------------------------------
x=sp.Symbol('x')
f=tanh(40*(x));
n_elements=5;
order=1;
psi=P(order);
domain=[-1, 1];

for e in range(0,n_elements+length(psi)-1):
    vertices[e]=domain[0] + e*(domain[1]-domain[0])/n_elements;


for e in range(0,n_elements):
    for s in range(0,order):
        cells[e, s]=e+s-1;


A = sp.zeros((n_elements, n_elements));
b_e = sp.zeros((n_elements, 1));


for e in range(0,n_elements):
    
    for s in range(0,order):
        omega[s]=vertices[cells[e, s]];

    
    for i in range(0,order):
        for j in range(0,order):
            A[cells[e, j], cells[e, j]]=A[cells[e, i], cells[e, j]]+calc_A(psi[i], psi[j], omega);

        b[cells[e, i]]=b[cells[e, i]]+calc_b(f, psi[i], omega);
    


c = A.LUsolve(b)
