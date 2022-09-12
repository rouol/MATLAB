%Нарисуй график количества итераций функции sinus в пределах от 0 до pi/2.
x = linspace(0, pi/2, 100);
[y, n] = sinus(x);
plot(x, n);
grid on;

