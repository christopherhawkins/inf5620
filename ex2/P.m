%% function to set up piecewise element functions
function [psi]=P(d)
d=d+1;
syms X;

domain=[-1 1];
for e=1:d
    Pa(e)=domain(1) + (e-1)*(domain(2)-domain(1))/(d-1);
end
for r=1:d
    psi(r)=0*X+1;

    for s=1:d
        if(r~=s)
           psi(r)=psi(r)*(X-Pa(s))/(Pa(r) - Pa(s));
        end
    end
end