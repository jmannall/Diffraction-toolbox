
close all
clear all
numPoints = 16;
plots = false;

x = [1, 0, 0];
test = x / vecnorm(x, 2, 2);

t_x = pi / 4;
t_y = 0;
t_z = pi / 4;

Rx = [1 0 0; 0 cos(t_x) -sin(t_x); 0 sin(t_x) cos(t_x)];
Ry = [cos(t_y) 0 sin(t_y); 0 1 0; -sin(t_y) 0 cos(t_y)];
Rz = [cos(t_z) -sin(t_z) 0; sin(t_z) cos(t_z) 0; 0 0 1];

mat = Rx * Rz;
test2 = (Rx*Rz*test')';

figure
plot3([0, test(1)], [0, test(2)], [0, test(3)])
hold on
plot3([0, test2(1)], [0, test2(2)], [0, test2(3)])
grid on
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])

points.Sun = Sunflower(numPoints, plots);
points.Fib = Fibinacci(numPoints, plots);

points.Tetrahedron = Tetrahedron(plots);
points.Octahedron = Octahedron(plots);
points.Cube = Cube(plots);
points.Icosahedron = Icosahedron(plots);
points.Dodecahedron = Dodecahedron(plots);

%% Combinations

plots = false;
numPoints = 32;

n = 0;
if numPoints < 4
elseif numPoints < 6
    points = Tetrahedron(plots, 0); % 4
elseif numPoints < 8
    points = Octahedron(plots); % 6
elseif numPoints < 12
    points = Cube(plots); % 8
elseif numPoints < 16
    points = Icosahedron(plots, 0); % 12
elseif numPoints < 20
    points = [Cube(false); Octahedron(plots)]; % 16
elseif numPoints < 24
    points = Dodecahedron(plots, 0); % 20
elseif numPoints < 32
    points = [Icosahedron(plots, 0); Icosahedron(plots, 1)]; % 24
else
    points = [Icosahedron(plots, 0); Dodecahedron(plots, 0)]; % 32
end

figure
plot3(0, 0, 0, 'o')
hold on
plot3(points(:,1), points(:,2), points(:,3), 'x')
title('Output')
grid on
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])

%% Sunflower

function points = Sunflower(numPoints, plot)
indices = 0:(numPoints - 1);

phi = acos(1 - 2*indices/numPoints);
theta = pi * (1 + 5^0.5) * indices;

x = cos(theta) .* sin(phi);
y = sin(theta) .* sin(phi);
z = cos(phi);
points = [x; y; z];

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(x, y, z, 'x')
    title('Sunflower')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

%% Fibinacci

function points = Fibinacci(numPoints, plot)
phi = pi * (sqrt(5) - 1);  % golden angle in radians

indices = 0:(numPoints - 1);
y = 1 - (indices / (numPoints - 1)) * 2; %  y goes from 1 to -1
radius = sqrt(1 - y .* y); % radius at y

theta = phi * indices; % golden angle increment

x = cos(theta) .* radius;
z = sin(theta) .* radius;
points = [x; y; z];

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(x, y, z, 'x')
    title('Fibinacci')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

%% Platonic solids

% Tetrahedron
function points = Tetrahedron(plot, idx)

points1 = [1, 1, 1;
    1, -1, -1;
    -1, 1, -1;
    -1, -1, 1];
points1 = points1 ./ vecnorm(points1, 2, 2);
points2 = -points1;
points2 = points2 ./ vecnorm(points2, 2, 2);

switch idx
    case 0
        points = points1;
    case 1
        points = points2;
    otherwise
        points = points1;
end

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(points1(:,1), points1(:,2), points1(:,3), 'x')
    plot3(points2(:,1), points2(:,2), points2(:,3), 'x')
    title('Tetrahedron')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

% Octahedron
function points = Octahedron(plot)

points = [1, 0, 0;
    0, 1, 0;
    0, 0, 1];
points = [points; -points];
points = points ./ vecnorm(points, 2, 2);

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(points(:,1), points(:,2), points(:,3), 'x')
    title('Octahedron')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

% Cube
function points = Cube(plot)

points = [1, 1, 1;
    1, -1, -1;
    -1, 1, -1;
    -1, -1, 1];
points = [points; -points];
points = points ./ vecnorm(points, 2, 2);

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(points(:,1), points(:,2), points(:,3), 'x')
    title('Cube')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

% Icosahedron
function points = Icosahedron(plot, idx)

phi = (1 + sqrt(5)) / 2;

points1 = [0, 1, phi;
    1, phi, 0;
    phi, 0, 1;
    0, -1, phi;
    -1, phi, 0;
    phi, 0, -1];
points1 = [points1; -points1];
points1 = points1 ./ vecnorm(points1, 2, 2);

points2 = [0, phi, 1;
    phi, 1, 0;
    1, 0, phi;
    0, phi, -1;
    phi, -1, 0;
    -1, 0, phi];
points2 = [points2; -points2];
points2 = points2 ./ vecnorm(points2, 2, 2);

switch idx
    case 0
        points = points2;
    case 1
        points = points1;
    otherwise
        points = points2;
end

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(points1(:,1), points1(:,2), points1(:,3), 'x')
    plot3(points2(:,1), points2(:,2), points2(:,3), 'x')
    title('Icosahedron')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

% Dodecahedron
function points = Dodecahedron(plot, idx)

phi = (1 + sqrt(5)) / 2;
invphi = 1 / phi;

points1 = [1, 1, 1;
    1, -1, -1;
    -1, 1, -1;
    -1, -1, 1;
    0, invphi, phi;
    invphi, phi, 0;
    phi, 0, invphi;
    0, -invphi, phi;
    -invphi, phi, 0;
    phi, 0, -invphi];
points1 = [points1; -points1];
points1 = points1 ./ vecnorm(points1, 2, 2);

points2 = [1, 1, 1;
    1, -1, -1;
    -1, 1, -1;
    -1, -1, 1;
    0, phi, invphi;
    phi, invphi, 0;
    invphi, 0, phi;
    0, phi, -invphi;
    phi, -invphi, 0;
    -invphi, 0, phi];
points2 = [points2; -points2];
points2 = points2 ./ vecnorm(points2, 2, 2);

switch idx
    case 0
        points = points1;
    case 1
        points = points2;
    otherwise
        points = points1;
end

if plot
    figure
    plot3(0, 0, 0, 'o')
    hold on
    plot3(points1(:,1), points1(:,2), points1(:,3), 'x')
    plot3(points2(:,1), points2(:,2), points2(:,3), 'x')
    title('Dodecahedron')
    grid on
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
end
end

