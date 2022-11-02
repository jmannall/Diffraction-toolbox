function [corners, planeCorners, planeRigid] = CreatePillar(i, j, gWidth, pWidth, height)
    step = (i - 1) * ((2 * j - 1) * (2 * pWidth + gWidth));
    corners = [pWidth step + pWidth 0
        -pWidth step + pWidth 0
        -pWidth step - pWidth 0
        pWidth step - pWidth 0
        pWidth step + pWidth height
        -pWidth step + pWidth height
        -pWidth step - pWidth height
        pWidth step - pWidth height];

    step = 8 * (2 * i - 3 + j);
    planeCorners = step + [1 4 3 2
        5 6 7 8
        1 2 6 5
        3 4 8 7
        2 3 7 6
        4 1 5 8];

    planeRigid = [0 0 1 1 1 1];
end   