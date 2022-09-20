%Numerical method of solving the wave problem.
%Differential equation: u_tt = c^2*u_xx
%Boundary conditions: u(0,t) = u(1,t) = 0
%Initial condition: u(x,0) = sin(2*pi*x), u_t(x,0) = 0
%Solution: u(x,t) = sin(2*pi*x)*cos(2*pi*t)


N = 1000;
h = 1/(N-1);
T = 2;
l = 1;
% tau <= h
tau = h;
M = ceil(T/tau);
x = linspace(0, l, N);
t = linspace(0, T, M);

c = 1;
a = c^2*tau^2/h^2;

u = zeros(N,M);
u(:,1) = sin(2*pi*x);
u_t_0 = zeros(N,1);
u(1,:) = 0;
u(N,:) = 0;

for n = 1:M-1
    if n == 1
        u(2:N-1,2) = u(2:N-1,1) + a/2*(u(3:N,1) - 2*u(2:N-1,1) + u(1:N-2,1)) + tau*u_t_0(2:N-1);
    else
        u(2:N-1,n+1) = 2*u(2:N-1,n) - u(2:N-1,n-1) + a*(u(3:N,n) - 2*u(2:N-1,n) + u(1:N-2,n));
    end
end

% analytical solution
u_an = sin(2*pi*x)'*cos(2*pi*t);

% graphs
figure(1)
colormap(jet);
mesh(x,t,u')
title('Численное решение')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')
figure(2)
colormap(jet);
mesh(x,t,u_an')
title('Аналитическое решение')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')

%Вывод ошибки
err = max(max(abs(u-u_an)));
disp(['Error = ',num2str(err)])
