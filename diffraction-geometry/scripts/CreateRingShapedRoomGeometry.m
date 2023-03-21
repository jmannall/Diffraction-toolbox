function [corners, planeCorners, source, receiver] = CreateRingShapedRoomGeometry(x, y, height, receiver)

    if nargin < 4
        receiver = [x(3) + 0.73 5.32 - 1.57 1.67];
    end
    corners = [0 0 0
        x(1) 0 0
        x(1) y(1) 0
        0 y(1) 0
        x(2) y(2) 0
        x(3) y(2) 0
        x(3) y(3) 0
        x(2) y(3) 0
        0 0 height
        x(1) 0 height
        x(1) y(1) height
        0 y(1) height
        x(2) y(2) height
        x(3) y(2) height
        x(3) y(3) height
        x(2) y(3) height];

      planeCorners = [1 2 3 4
          9 12 11 10
          1 9 10 2
          2 10 11 3
          3 11 12 4
          4 12 9 1
          5 6 14 13
          6 7 15 14
          7 8 16 15
          8 5 13 16];

    source = [x(2) - 0.8 y(1) - 18.39 - 1.84 1.46];
end