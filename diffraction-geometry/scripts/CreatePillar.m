function [corners, planeCorners, planeRigid] = CreatePillar(i, j, width, height)
    step = 3 * (2 * j - 1) * (i - 1) * width;
    corners = [step + width width 0
        step + width -width 0
        step - width -width 0
        step - width width 0
        step + width width height
        step + width -width height
        step - width -width height
        step - width width height];

    step = 8 * (2 * i - 3 + j);
    planeCorners = step + [1 2 3 4
        5 8 7 6
        1 5 6 2
        3 7 8 4
        2 6 7 3
        4 8 5 1];

    planeRigid = [0 0 1 1 1 1];
end   