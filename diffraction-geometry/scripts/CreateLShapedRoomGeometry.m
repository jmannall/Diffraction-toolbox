function [corners, planeCorners, source, receiver] = CreateLShapedRoomGeometry(x, y, height, receiver)

    if nargin < 4
        receiver = [sum(x) sum(y) height] / 2;
    end
    corners = [0 0 0
        x(1) 0 0
        x(1) y(1) 0
        x(2) y(1) 0
        x(2) y(2) 0
        0 y(2) 0
        0 0 height
        x(1) 0 height
        x(1) y(1) height
        x(2) y(1) height
        x(2) y(2) height
        0 y(2) height];

      planeCorners = [1 2 3 4 5 6
          7 12 11 10 9 8
          1 7 8 2 0 0
          3 9 10 4 0 0
          2 8 9 3 0 0
          4 10 11 5 0 0
          5 11 12 6 0 0
          6 12 7 1 0 0];

%             planeCorners = [1 13 5 6
%           7 12 11 14
%           1 7 8 2
%           3 9 10 4
%           2 8 9 3
%           4 10 11 5
%           5 11 12 6
%           6 12 7 1];

    source = [x(1) y(1) height] / 2;
end