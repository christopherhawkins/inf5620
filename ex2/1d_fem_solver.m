clear;

%% define system to solve
syms x; % symbol to be used in equations
f=tanh(40*(x)); % function to solve
n_elements=1; % number of elements
order=4;
psi=P(order); % type of piecewise functions
domain=[-1 1]; % domain of the system
name='Exercise 15 P4 1';

%% set up vertices and cells
for e=1:n_elements*(length(psi)-1)+1
    vertices(e)=domain(1) + (e-1)*(domain(2)-domain(1))/(n_elements*(length(psi)-1));
end

for e=1:n_elements
    for s=1:length(psi)
        cells(e, s)=((e-1)*(length(psi)-1))+s;
    end
end

%% set up matrices A and b
A=zeros(length(vertices), length(vertices));
b=0*(1:length(vertices));

%% solve system
% loop over all elements
for e=1:length(cells(:, 1))
    
% set up local domain
    for s=1:length(cells(e, :))
        omega(s)=vertices(cells(e, s));
    end
    
% calculate A and b
    for i=1:length(psi)
        for j=1:length(psi)
            A(cells(e, i), cells(e, j))=A(cells(e, i), cells(e, j))+calc_A(psi(i), psi(j), omega);
        end
        b(cells(e, i))=b(cells(e, i))+calc_b(f, psi(i), omega);
    end
    
end

%% calculate coefficients
m(:, 1)=b;
c=A\m;

%% plot result
x_plot=domain(1):0.01:domain(2);
h=figure();
hold on;

p1=ezplot(f,domain);
set(p1,'Color','k');

syms X;

for e=1:n_elements
    element_plot(e)=X*0.0001;
    for s=1:length(psi)
        omega(1)=vertices(cells(e, 1));
        omega(2)=vertices(cells(e, length(psi)));
        psi_local(s)=subs(psi(s), 2*(X-omega(1))/(omega(2)-omega(1)) -1 );
        element_plot(e)=element_plot(e)+c(cells(e, s))*psi_local(s);
    end
    ezplot(element_plot(e), [vertices(cells(e, 1)) vertices(cells(e, length(psi)))]);
end

axis([-1 1 -2 2])
title(['P' num2str(order) ' functions, fitted to f=tanh(40*(x)) using ' num2str(n_elements) ' elements'])
xlabel('x')
legend('tanh(40*(x))', 'FEM')
print(h,'-dpdf', [name '.pdf'])

