function [room, source, receiver] = CreateRingShapedRoomGeometry(x, y, z)

    source = [(x(3) + x(1)) / 2, (y(3) + y(2)) / 2, min(1.6, z - 0.1)];
    receiver = [x(2) / 2, (y(3) + y(2)) / 2, min(1.6, z - 0.1)];

    room.corners = [0 0 0
        x(1) 0 0
        x(1) y(1) 0
        0 y(1) 0
        x(2) y(2) 0
        x(3) y(2) 0
        x(3) y(3) 0
        x(2) y(3) 0
        0 0 z
        x(1) 0 z
        x(1) y(1) z
        0 y(1) z
        x(2) y(2) z
        x(3) y(2) z
        x(3) y(3) z
        x(2) y(3) z];

    room.planeCorners = [1 2 3 4
        9 12 11 10
        1 9 10 2
        2 10 11 3
        3 11 12 4
        4 12 9 1
        5 6 14 13
        6 7 15 14
        7 8 16 15
        8 5 13 16];

    room.edgeCorners = [5 13
        6 14
        7 15
        8 16];

    room.edgeCanSeeEdge = [false true false true
        true false true false
        false true false true
        true false true false];

    room.thetaW = 270 * ones(4, 1);
end