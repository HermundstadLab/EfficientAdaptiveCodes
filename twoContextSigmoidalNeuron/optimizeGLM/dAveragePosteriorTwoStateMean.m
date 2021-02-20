function dAv = dAveragePosteriorTwoStateMean(x, y, pLowPrev, muL, muH, h)
%Compute a distortion measure as a squared difference of average posteriors
%for two-state mean environment with a fixed variance equal to 1
%Input:
%   x, y - quantized variable and output variable
%   pLowPrev - prior on low variance
%   muL (optional) - low mean - default 1
%   muH (optional) - high mean - default 2
%   h (optional) - hazard rate - default 1e-2

if nargin < 4
    muL = 1;
    muH = 2;
    h = 1e-2;
end

%compute average posterior for x
pLowX = (((1-h)*pLowPrev + h*(1-pLowPrev)) / sqrt(2 * pi)) .* exp(-(x - muL).^2 ./ 2);
Zx = pLowX + (((h*pLowPrev + (1-h) * (1-pLowPrev)) / sqrt(2 * pi)) * exp(-(x - muH).^2 ./ 2));
pLowX = 1 ./ Zx .* pLowX;
pHighX = 1 - pLowX;
avgX = pLowX .* muL + pHighX .* muH;

%compute average posterior for y
pLowY = (((1-h)*pLowPrev + h*(1-pLowPrev)) / sqrt(2 * pi)) * exp(-(y - muL).^2 ./ 2);
Zy = pLowY + (((h*pLowPrev + (1-h) * (1-pLowPrev)) / sqrt(2 * pi)) * exp(-(y - muH).^2./2));
pLowY = 1 ./ Zy .* pLowY;
pHighY = 1 - pLowY;
avgY = pLowY .* muL + pHighY .* muH;

dAv = (avgX - avgY).^2;
