% heat equation for rode with two parts of different material and two ends of different temperature
% 1D heat equation
% du/dt = alpha * d2u/dx2, alpha = k/(rho*c)
clc
clear all

L1 = 0.01;
L2 = 0.09;
L = L1 + L2;

k1 = 400;
k2 = 237;
rho1 = 8940;
rho2 = 2700;
c1 = 385;
c2 = 920;
alpha1 = k1/(rho1*c1);
alpha2 = k2/(rho2*c2);

T0 = 20;
TA = 120;
TC = 35;

% 1D heat equation
% du/dt = alpha * d2u/dx2, alpha = k/(rho*c)
% u(x,0) = T0
% u(0,t) = TA
% u(L,t) = TC

% analytical solution
TB = ((k1/L1)*TA + (k2/L2)*TC) / ((k1/L1) + (k2/L2));
u1 = @(x) (TB-TA)*x/L1 + TA;
u2 = @(x) (TC-TB)*(x-L1)/L2 + TB;
u_true = @(x) (x<=L1).*u1(x) + (x>L1).*u2(x);

% plot analytical solution
x = linspace(0,L,100);
plot(x,u_true(x),'r-','LineWidth',2);
xlabel('x');
ylabel('u');
hold on;

% numerical solution
M = 200;
N = 1800;
K = 200;
h = (L1 + L2) / (M+N);
T = 30;
tau = T / K;

x = linspace(0, L1 + L2, M + N + 1);
t = linspace(0, T, K + 1);

u = zeros(M + N + 1, K + 1);
u(:, 1) = T0;
u(1, :) = TA;
u(M + N + 1, :) = TC;

% time loop
for k = 1:K
    % plot numerical solution
    plot(x, u(:, k), 'b-');
    pause(0.1);
    % stop if solution is close to analytical solution with tolerance 0.01
    if (mean(abs(u(:, k) - u_true(x)')) < 0.5)
        time = k * tau
        break;
    end

    A = zeros(M + N - 1, M + N - 1);
    b = zeros(M + N - 1, 1);
    % x = 0
    A(1, 1) = 1;
    b(1) = TA;
    % x < L1
    for m = 2:M
        A(m, m - 1) = alpha1 * tau / h^2;
        A(m, m) = -2 * alpha1 * tau / h^2 - 1;
        A(m, m + 1) = alpha1 * tau / h^2;
        b(m) = -u(m, k);
    end
    % x = L1
    A(M, M - 1) = -k1;
    A(M, M) = k1 + k2;
    A(M, M + 1) = -k2;
    b(M) = 0;
    % x > L1
    for m = M + 1:M + N - 2
        A(m, m - 1) = alpha2 * tau / h^2;
        A(m, m) = -2 * alpha2 * tau / h^2 - 1;
        A(m, m + 1) = alpha2 * tau / h^2;
        b(m) = -u(m, k);
    end
    % x = L
    A(M + N - 1, M + N - 1) = 1;
    b(M + N - 1) = TC;
    u(2:M + N, k + 1) = A \ b;
end

% T = 100
% N = 200;
% dx = L/N;
% dt = 0.01;

% initial condition
% u = zeros(N+1,1);
% u(1) = TA;
% u(2:N) = T0;
% u(N+1) = TC;

% % make function of alpha using heaviside function
% syms x
% % alpha = @(x) alpha2 + (alpha1 - alpha2) * heaviside(x-L1);
% alpha = alpha2 + (alpha1 - alpha2) * heaviside(x-L1);
% sym_alpha = symfun(alpha, x);
% deriv_alpha = diff(sym_alpha, x);
% % time loop
% for t = 1:2
%     % space loop
%     for i = 2:N
%         i*dx
%         u(i) = u(i) + double(sym_alpha(i*dx))*dt/dx^2 * (u(i+1) - 2*u(i) + u(i-1)) + double(deriv_alpha(i*dx))*dt/dx^2 * (u(i+1) - u(i-1));
%     end
%     % plot numerical solution
%     plot(dx*(0:N),u,'b-');
%     t
%     pause(0.01);
% end
% double(sym_alpha(0))
% double(deriv_alpha(0.01))