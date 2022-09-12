x = linspace(25, 50, 1001);

y1 = sinus(x);
y2 = sinus_plus(x);

plot(x, y1, x, y2);

xlabel('x');
ylabel('y');
legend('sinus', '\_sinus\_plus');
grid on;
