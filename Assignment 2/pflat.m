function [y] = pflat(x)
    y = [x(1:end-1,:) ./ x(end,:); ones(1, size(x,2))];
end

