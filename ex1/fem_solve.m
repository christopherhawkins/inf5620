function[c]=fem_solve(f, psi, omega)

N=length(psi)-1;

%% Set up matricies A and b (A=cb)
b = zeros(N, 1); % matrix of RHS (i.e. b)
A = zeros(N, N);   % local K matrix

%% Calculate variables for matrix A
for i=0:N
    for j=0:N
        A(i+1, j+1)=int(psi(i+1)*psi(j+1), omega(1), omega(2));
    end
end

%% Calculate variables for matrix b
for i=0:N
    b(i+1)=int(f*psi(i+1), omega(1), omega(2));
end

%% Divide matricies to finc coefficients c
c=A\b;