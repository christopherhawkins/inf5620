clear;

%% Define variables to solve equation x using power law series x^i to solve f
syms x; % symbol to be used in equations
omega=[0 2*pi]; % size of domain from omega(1) -> omega(2)

%% Define equation to solve
s=20;
f=tanh(s*(x-pi));

% loop over a number of N values
for N=3:10;
    %% Set up series psi
    clear psi;
    for i=0:N;
        psi(i+1)=sin((2*i+1)*x);
    end

    %% Solve
    c=fem_solve(f, psi, omega);

    %% Create plotting variables for the solution
    xp=(omega(2)-omega(1))*(1:1000)/(1000);
    for i=1:1000;
        solution(i)=0;
        for j=0:N;
            solution(i)=solution(i)+c(j+1)*sin((2*j+1)*xp(i));
        end
    end

    %% Plot result
    clf;
    hold on;
    plot(xp, tanh(s*(xp-pi)), 'k');
    plot(xp, solution, '--', 'color', [0.5 0.5 0.5]);
    xlabel('x');
    ylabel('f(x)')
    title(['N=' num2str(N)])
    legend('f(x)', 'solution')
    drawnow;
end
