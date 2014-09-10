clear;

%% Define variables to solve equation x using power law series x^i to solve f
syms x; % symbol to be used in equations
N=[2 4 6]; % number of terms in power series
omega=[0 8]; % size of domain from omega(1) -> omega(2)

%% Define equation to solve
f=exp(-x);

%% Set c and loop over all power law lengths N
c=zeros(7, 3);
for s=1:length(N)
    %% Set up power law series
    clear psi;
    for i=0:N(s);
        psi(i+1)=x^i;
    end

    %% Solve
    ci=fem_solve(f, psi, omega);
    c(1:length(ci), s)=ci(:);
end

%% Create plots for each selected N
for s=1:length(N)
    %% Create plotting variables for the solution
    x=(omega(2)-omega(1))*(1:1000)/(1000);
    for i=1:1000;
        solution(i)=0;
        for j=0:N(s);
            solution(i)=solution(i)+c(j+1, s)*x(i)^j;
        end
    end

    %% Plot result
    figure();
    hold on;
    plot(x, exp(-x), 'k');
    plot(x, solution, '--', 'color', [0.5 0.5 0.5]);
    xlabel('x');
    ylabel('f(x)')
    title(['N=' num2str(N(s))])
    legend('f(x)', 'solution')

end
