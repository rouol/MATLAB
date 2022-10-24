% Построение анимационного графика функции Вейерштрасса
% x - переменная, t - время
% w - входной параметр < 1
% x лежит в диапазоне [-10,10]
% t лежит в диапазоне [0,2*pi]

function AnimatePlot(n)
    if n < 1
        error('n must be more than 1');
    end
    % Построение анимации в каждый момент времени t c сохранением в gif-файл
    % Параметры анимации
    t = 0:0.01:2*pi;
    x = -10:0.1:10;
    
    % Параметры анимации
    fps = 30; % Количество кадров в секунду
    delay = 1/fps; % Задержка между кадрами
    filename = 'Weierstrass.gif'; % Имя файла для сохранения анимации

    a = 1/2;
    b = 3;
    
    % Построение анимации
    for a = 0.05:.05:1
        y = zeros(length(x),5);
        [row, col] = size(y);
        Y = zeros(size(x));
        for i = 1:length(x)
            for j = 1:col
                y(i,j) = (a^(j-1))*cos((b^(j-1))*pi*x(i));
            end
            Y(i) = sum(y(i,:));
        end
        plot(x,Y);
        axis([-10 10 -1.5 1.5]);
        title(['a = ', num2str(a)]);
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if a == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',delay);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay);
        end
        % pause(0.1);
    end
end