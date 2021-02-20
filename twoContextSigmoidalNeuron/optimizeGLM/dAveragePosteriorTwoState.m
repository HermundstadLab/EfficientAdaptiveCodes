function dAv = dAveragePosteriorTwoState(x, y, pLowPrev, sigmaL, sigmaH, h)
%Compute a distortion measure as a squared difference of average posteriors
%for two-state variance environment
%Input:
%   x, y - quantized variable and output variable
%   pLowPrev - prior on low variance
%   sigmaL (optional) - low standard deviation - default 1
%   sigmaH (optional) - high standard deviation - default 2
%   h (optional) - hazard rate - default 1e-2

if nargin < 4
    sigmaL = 1;
    sigmaH = 2;
    h = 1e-2;
end

%compute average posterior for x
pLowX = (((1-h)*pLowPrev + h*(1-pLowPrev)) / sqrt(2 * pi * sigmaL^2)) .* exp(-x.^2./(2*sigmaL^2));
Zx = pLowX + (((h*pLowPrev + (1-h) * (1-pLowPrev)) / sqrt(2 * pi * sigmaH^2)) * exp(-x.^2./(2*sigmaH^2)));
pLowX = 1 ./ Zx .* pLowX;
pHighX = 1 - pLowX;
avgX = pLowX .* sigmaL + pHighX .* sigmaH;

%compute average posterior for y
pLowY = (((1-h)*pLowPrev + h*(1-pLowPrev)) / sqrt(2 * pi * sigmaL^2)) * exp(-y.^2./(2*sigmaL^2));
Zy = pLowY + (((h*pLowPrev + (1-h) * (1-pLowPrev)) / sqrt(2 * pi * sigmaH^2)) * exp(-y.^2./(2*sigmaH^2)));
pLowY = 1 ./ Zy .* pLowY;
pHighY = 1 - pLowY;
avgY = pLowY .* sigmaL + pHighY .* sigmaH;

dAv = (avgX - avgY).^2;
