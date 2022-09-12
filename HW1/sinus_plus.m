%Значение синуса вычисляются от 0 до 2pi по формуле ряда Тейлора, а в других диапазонах используется периодичность синуса.
function y = sinus_plus(x)
    y = zeros(size(x));
    for i = 1:length(x)
        if x(i) >= 0 && x(i) <= 2*pi
            y(i) = sinus(x(i));
        else
            y(i) = sinus_plus(x(i) - 2*pi*floor(x(i)/(2*pi)));
        end
    end
end