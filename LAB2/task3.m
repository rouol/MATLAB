% %x = linspace(-2,2);
% x = -2:.004:2;
% %a = 1/4;
% %b = 51;
% a = 1/2;
% b = 3;
% y = zeros(length(x),5);
% [row, col] = size(y);
% Y = zeros(size(x));

% for i = 1:length(x)
%     for j = 1:col
%         y(i,j) = (a^(j-1))*cos((b^(j-1))*pi*x(i));
%     end
%     Y(i) = sum(y(i,:));
% end
% plot(x,Y)

AnimatePlot(10)