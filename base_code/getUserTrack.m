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

function [track,velocity_dir] = getUserTrack(TrackType, InitPos, distance, d_update, direction, varargin)

% 'direction' is the azimuthal angle between the north and the track clockwisely 
% Input for 'Hexagon'-type track:
%       'sidelength': the length of the side of a hexagon
%       'Orientation': '0' means counter-clockwise; '1' means clockwise
%
% Output: The snapshot positions on the track: 1 m separation distance
%         between two snapshots.

if ~strcmp(TrackType, 'Linear') && ~strcmp(TrackType, 'Hexagon')
    error('The input track type is not supported yet.');
end


tr_length = ceil(distance);
% d_snap = 1;
no_dps = tr_length/d_update;
switch TrackType
    case 'Linear'
%         if distance > 100 || distance < 0
        if distance < 0    
            error('The track length is out of range.')
        end
        
        RelativePos = [cos(deg2rad(direction));sin(deg2rad(direction));0]*(0:(no_dps-1))*d_update;
        track = InitPos + RelativePos;
        velocity_dir = ones(1,no_dps)*direction;
        
        
    case 'Hexagon'
        var_names = {'side_length', 'orientation'};
        var_defaults = [5,1];
        
        for n = 1:2
            if numel( varargin ) >= n && ~isempty( varargin{n} )
                if ~(isnumeric( varargin{n} ) && varargin{n}>=0 &&...
                        isreal( varargin{n} ) && all(size(varargin{n}) == [1 1]))
                    error([var_names{n},' has wrong format']);
                end
                
                eval([ var_names{n},'=',num2str(num2str(varargin{n})),';'  ]);
            else
                eval([ var_names{n},'=',num2str(var_defaults(n)),';'  ]);
            end
        end
        side_snap = side_length/d_update; % # of snapshots on each side

        if orientation == 1
            theta = fliplr(0:pi/3:2*pi)+2*pi/3;
            velocity_dir = [ones(1,side_snap)*direction,...
            ones(1,side_snap)*mod((direction-60),360),...
            ones(1,side_snap)*mod((direction-120),360),...
            ones(1,side_snap)*mod((direction-180),360),...
            ones(1,side_snap)*mod((direction-240),360),...
            ones(1,side_snap)*mod((direction-300),360)];
        elseif orientation == 0
            theta = (0:pi/3:2*pi)-2*pi/3;
            velocity_dir = [ones(1,side_snap)*direction,...
            ones(1,side_snap)*mod((direction+60),360),...
            ones(1,side_snap)*mod((direction+120),360),...
            ones(1,side_snap)*mod((direction+180),360),...
            ones(1,side_snap)*mod((direction+240),360),...
            ones(1,side_snap)*mod((direction+300),360)];
        else
            error('Wrong input orientation');
        end
        RelativePos = zeros(3,6*side_snap);
        % rotation angle east is 0
        rot_matrix = [cos(deg2rad(direction)) -sin(deg2rad(direction));...
            sin(deg2rad(direction)) cos(deg2rad(direction))];
        x_hex = cos(theta);
        y_hex = sin(theta);
%         coor_hex = [x_hex;y_hex]*side_length;
        
        coor_hex = [x_hex;y_hex]*side_length;
        
        coor_rot = rot_matrix*coor_hex; % rotate the hexagon based on direction
        coor_trans = coor_rot - coor_rot(:,1); % Move the start point to the origin

        for i = 1:6
            st = coor_trans(:,i);
            en = coor_trans(:,i+1);
            sep = (en-st)/side_snap*(0:(side_snap-1));
            pos = st+sep;
            RelativePos(1:2,(i-1)*side_length/d_update+1:i*side_length/d_update)...
                = pos; 
        end
        circles = ceil(no_dps/6/side_snap);
        RelativePos = repmat(RelativePos,1,circles);
        track = InitPos+RelativePos(:,1:no_dps)*d_update;

        velocity_dir = repmat(velocity_dir,1,circles);
        velocity_dir = velocity_dir(1:no_dps);
end
        
        
            