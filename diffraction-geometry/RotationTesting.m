
% Rotates clockwise. Add a translation part and cna create function to
% rotate sourecs and receivers around edges
x = 0;
y = 0;
z = 1;

u = [x y z];
u = u / norm(u);

theta = 90;
c = cosd(theta);
s = sind(theta);
C = 1 - c;

d = [c -z * s y * s
    z * s c -x * s
    -y * s x * s c];
Q = C * (u' * u) + d;

point = [1 0 1];
rotPoint = point * Q;