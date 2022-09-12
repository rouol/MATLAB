%график sinus и sin на одном графике в пределах от 25 до 50 с шагом 0.01

x = 25:0.01:50;
y = sinus(x);

plot(x,y)
hold on

plot(x,sin(x))
hold off

grid on

xlabel('x')
ylabel('y')

legend('sinus','sin')

