function [y, n] = sinus(x)
    y = zeros(1, length(x));
    n = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = 0;
        n(i) = 0;
        while 1
            n(i) = n(i) + 1;
            y(i) = y(i) + ((-1)^(n(i) - 1) * x(i)^(2 * n(i) - 1)) / factorial(2 * n(i) - 1);
            if abs(((-1)^(n(i) - 1) * x(i)^(2 * n(i) - 1)) / factorial(2 * n(i) - 1)) == 0
                break
            end
        end
    end
end



