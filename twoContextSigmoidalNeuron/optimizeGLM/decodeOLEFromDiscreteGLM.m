function [xHat, params] = decodeOLEFromDiscreteGLM(y, x)
%Optimal Linear Estimator (OLE) from GLM output

if numel(unique(y))<2
    xHat = zeros(size(y));
    params = [0, 0];
else
    params = polyfit(y, x, 1);
    xHat = y * params(1) + params(2);
end


