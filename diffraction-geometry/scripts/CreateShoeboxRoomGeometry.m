function room = CreateShoeboxRoomGeometry(x, y, z)

    room.receiver = [x / 2, y / 4, min(1.6, z - 0.1)];
    room.source = [x / 2, 3 * y / 4, min(1.6, z - 0.1)];

    room.corners = [0 0 0
        x(1) 0 0
        x(1) y(1) 0
        0 y(1) 0
        0 0 z
        x(1) 0 z
        x(1) y(1) z
        0 y(1) z];

    room.planeCorners = [1 2 3 4
        5 8 7 6
        1 5 6 2
        2 6 7 3
        3 7 8 4
        4 8 5 1];
end