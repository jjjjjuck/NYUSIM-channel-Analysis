%%% NYUSIM - User License %%%

% Copyright (c) 2016-2019 New York University and NYU WIRELESS

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the “Software”),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software. Users shall cite 
% NYU WIRELESS publications regarding this work.

% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUTWARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
% OTHER LIABILITY, WHETHER INANACTION OF CONTRACT TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.

function LosMap = getLosMap(area,d_co,h_BS,h_MS,sceType)

% Generate the map of spatially correlated LOS/NLOS condiiton

% Generate grid
d_px = 1;
N = area/d_px + 1;
delta = d_co/d_px;

[x,y] = meshgrid(-area/2:d_px:area/2);
Pr_LOS = zeros(N,N);

d_2D = sqrt(x.^2+y.^2)+eps;
d1 = 22;
d2 = 100;

% Probability of LOS/NLOS condition (depends on model and scenario)
for i = 1:N
    for j = 1:N
        switch sceType
            % NYU squared model for 5GCM
            case 'UMi'
                Pr_LOS(i,j) = (min(d1/d_2D(i,j),1)*(1-exp(-d_2D(i,j)/d2))+...
            exp(-d_2D(i,j)/d2))^2;
            % NYU squared model for 5GCM
            case 'UMa'
                if h_MS < 13 || h_BS <= 18
                    C = 0;
                elseif h_MS > 13 && h_MS < 23
                    C = ((h_MS-13)/10)^1.5*1.25e-6*d_2D(i,j)^3*exp(-d_2D(i,j)/150);
                else
                    disp('Height of base station or mobile terminal is out of range.');
                end
                Pr_LOS(i,j) = ((min(d1/d_2D(i,j),1)*(1-exp(-d_2D(i,j)/d2))+...
            exp(-d_2D(i,j)/d2))*(1+C))^2;
            % 3GPP 38.901 Release 14
            case 'RMa'
                if d_2D(i,j) <= 10
                    Pr_LOS(i,j) = 1;
                elseif d_2D(i,j) > 10
                    Pr_LOS(i,j) = exp(-(d_2D(i,j)-10)/1000);
                end
        end
    end
end

% The filter coefficient is considered as 0 when the distance is beyond
% 4*d_co.
M = 8*delta+1;
h = zeros(M);
InitMap = randn(N+M-1);

% Generate the filter
for i = 1:M
    for j = 1:M
        h(i,j) = exp(-sqrt(((M+1)/2-i).^2+((M+1)/2-j).^2)/d_co);
    end
end
% h((M+1)/2,(M+1)/2) = 0;
CorrMapPad = conv2(InitMap,h,'same');
CorrQ = CorrMapPad((M+1)/2:(M+1)/2+N-1,(M+1)/2:(M+1)/2+N-1);
CorrK = 1/2*(1+erf(CorrQ/sqrt(2)));
LosMap = double(CorrK < Pr_LOS);

end


