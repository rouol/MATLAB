%Numerical method for solving the thermal conductivity problem.
%Differential equation: u_t = u_xx
%Boundary conditions: u(0,t) = u(1,t) = 0
%Initial condition: u(x,0) = sin(2*pi*x)
%Solution: u(x,t) = sin(2*pi*x)*exp(-4*pi^2*t)

a = 0; b = 1;
N = 100;
h = (b-a)/N;

lambda = 0.5;

T = 0.1;
tau = lambda*h^2;
M = round(T/tau);

u = zeros(N,M);
x = linspace(a,b,N);
t = linspace(0,T,M);

% initial condition
u(:,1) = sin(2*pi*x);

% boundary conditions
u(1,:) = 0;
u(N,:) = 0;

% numerical solution
for n = 1:M-1
    u(2:N-1,n+1) = u(2:N-1,n) + tau/h^2*(u(3:N,n)-2*u(2:N-1,n)+u(1:N-2,n));
end

% numerical solution graph
% change color map
colormap(jet);
mesh(x,t,u');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Численное решение');

% analytical solution
u_an = sin(2*pi*x)'*exp(-4*pi^2*t);

% analytical solution graph
figure;
% change color map
colormap(jet);
mesh(x,t,u_an');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
title('Аналитическое решение');

% error output
disp('Error:');
disp(max(max(abs(u_an-u))));

% % 60 fps gif animation
% figure;
% for n = 1:M
%     plot(x,u(:,n));
%     axis([a b -1 1]);
%     xlabel('x');
%     ylabel('u(x,t)');
%     title(['Numerical solution, t = ',num2str(t(n))]);
%     drawnow;
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if n == 1
%         imwrite(imind,cm,'animation.gif','gif','Loopcount',inf, 'DelayTime',1/60);
%     else
%         imwrite(imind,cm,'animation.gif','gif','WriteMode','append', 'DelayTime',1/60);
%     end
% end
