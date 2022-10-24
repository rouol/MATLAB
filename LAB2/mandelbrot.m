function mandelbrot(R, n)

    x0 = -2;   x1 = 1;
    y0 = -1.5; y1 = 1.5;
    
    [x,y] = meshgrid(linspace(x0, x1, R), linspace(y0, y1, R));
    
    c = x + 1i * y;
    z = zeros(size(c));
    k = zeros(size(c));
    
    for ii = 1:n
        z   = z.^2 + c;
        k(abs(z) > 2 & k == 0) = n - ii;
    end
    
    figure,
    imagesc(k),
    colormap hsv
    axis square