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

function init_pos = getInitPos(subpath_AODs, subpath_AOAs,TRDistance, envType, h_MS)

% Obtain the initial position of the user terminal

phi_AOD = subpath_AODs.c1.SP1(1);
theta_ZOD = subpath_AODs.c1.SP1(2);

if strcmp(envType,'LOS')
    init_pos = [TRDistance*cos(deg2rad(theta_ZOD))*sin(deg2rad(phi_AOD));
    TRDistance*cos(deg2rad(theta_ZOD))*cos(deg2rad(phi_AOD));
    h_MS];
elseif strcmp(envType,'NLOS')
    phi_AOA = subpath_AOAs.c1.SP1(1);
    phi_AOA = mod(phi_AOA+180,360);
    dz = 1/TRDistance;
    ang1 = min(deg2rad(phi_AOA),deg2rad(phi_AOD));
    ang2 = max(deg2rad(phi_AOA),deg2rad(phi_AOD));
    if ang2-ang1 <= pi
        ang = ang1:dz:ang2;
    else
        ang = ang2:dz:(ang1+2*pi);
    end
    pt = randi(length(ang)-2);
    ang_TR = ang(pt+1);
    init_pos = [TRDistance*cos(deg2rad(theta_ZOD))*sin(ang_TR);
    TRDistance*cos(deg2rad(theta_ZOD))*cos(ang_TR);
    h_MS];
end