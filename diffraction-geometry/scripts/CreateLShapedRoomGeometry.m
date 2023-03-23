function [room, source, receiver] = CreateLShapedRoomGeometry(x, y, z)

    source = [(x(2) + x(1)) / 2, (y(2) + y(1)) / 2, min(1.6, z - 0.1)];
    receiver = [x(1) / 2, y(1) / 2, min(1.6, z - 0.1)];

    room.corners = [0 0 0
        x(1) 0 0
        x(1) y(1) 0
        x(2) y(1) 0
        x(2) y(2) 0
        0 y(2) 0
        0 0 z
        x(1) 0 z
        x(1) y(1) z
        x(2) y(1) z
        x(2) y(2) z
        0 y(2) z];

    room.planeCorners = [1 2 3 4 5 6
        7 12 11 10 9 8
        1 7 8 2 0 0
        3 9 10 4 0 0
        2 8 9 3 0 0
        4 10 11 5 0 0
        5 11 12 6 0 0
        6 12 7 1 0 0];

    room.edgeCorners = [3 9];

    room.thetaW = 270;
end