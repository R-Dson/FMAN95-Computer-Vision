
%% load 1)
load('CompEx1.mat');

% projection

% The 2d projection on (x1, x2, 1) is given by 
% (x1, x2) = (x'1/x'3, x'2/x'3) with the extra value 1
%pflat = @(x) [x(1:end-1,:) ./ x(end,:); ones(1, size(x,2))];

%% Computer Exercise 1
% lambda variable
syms lambda;

x1 = [4; -2; 2];
x2 = [3; -2; -1];
x3 = lambda .* [4; -2; 2];

% Just a test that plots a 3d cube on the 2d plane
% d = [-1 -1 -1 -1 1 1 1 1; -1 -1 1 1 -1 -1 1 1; 2 4 2 4 2 4 2 4];
% t = pflat(d);
% plot3(t(1,:),t(2,:),t(3,:), '.')


v1 = pflat(x1)
v2 = pflat(x2)
v3 = pflat(x3)

plot(v1, '*')
hold on;
plot(v2, '*')
plot(v3, '*')

%% 
v4 = pflat(x2D);
v5 = pflat(x3D);
plot3(v4(1,:),v4(2,:),v4(3,:), '.')
hold on;
plot3(v5(1,:),v5(2,:),v5(3,:), '.')

%% Computer Exercise 2
y = imread('compEx2.JPG');
figure
colormap gray
imagesc(y)
hold on;

load('CompEx2.mat');

% plotting the points
plot(p1(1,:), p1(2,:), 'x', 'color','r');
plot(p2(1,:), p2(2,:), 'x', 'color','r');
plot(p3(1,:), p3(2,:), 'x', 'color','r');

% using null to solve for the line and transpose (l^T*x=0)
lc1 = null(p1');
lc2 = null(p2');
lc3 = null(p3');

rital(lc1);
rital(lc2);
rital(lc3);
%
% I tried to do this at first but it was off by a few pixels
% constants = [lc1(end ); lc2(end); lc3(end )];
% ab = [lc1(1:end-1 )'; lc2(1:end-1 )'; lc3(1:end-1 )'];
% inter = abs(ab\constants)
% plot(inter)

% Solving the equation l_1^T * x=0 and l_2^T * x=0 in R3
ns = null([lc2 lc3]');
% To the plane P2
P = pflat(ns);
plot(P(1,:), P(2,:), 'x', 'color','r');
title('');
% Calculating the distance 
d = abs(lc1' * P) / sqrt(lc1(1)^2 + lc1(2)^2)

%% Computer Exercise 3

load('CompEx3.mat');

H1 = [sqrt(3) -1 1; 1 sqrt(3) 1; 0 0 2];
H2 = [1 -1 1; 1 1 0; 0 0 1];
H3 = [1 1 0; 0 2 0; 0 0 1];
H4 = [sqrt(3) -1 1; 1 sqrt(3) 1; 1/4 1/2 2];

% adding a 1 at the end of each point
sp = [startpoints; ones(1, size(startpoints, 2))];
ep = [endpoints; ones(1, size(startpoints, 2))];

% multiplying and mapping to the plane
pflats1 = pflat(H1*sp);
pflate1 = pflat(H1*ep);

pflats2 = pflat(H2*sp);
pflate2 = pflat(H2*ep);

pflats3 = pflat(H3*sp);
pflate3 = pflat(H3*ep);

pflats4 = pflat(H4*sp);
pflate4 = pflat(H4*ep);

% original
figure;
pflats = pflat(sp);
pflate = pflat(ep);
plot ([ pflats(1 ,:); pflate(1 ,:)] , ...
[pflats(2 ,:); pflate(2 ,:)], 'b - ');

%% plotting the figures 
figure;
t = tiledlayout(2,2);
nexttile
plot ([ pflats1(1 ,:); pflate1(1 ,:)] , ...
[pflats1(2 ,:); pflate1(2 ,:)], 'b - ');
title('H1')
axis equal

nexttile
plot ([ pflats2(1 ,:); pflate2(1 ,:)] , ...
[pflats2(2 ,:); pflate2(2 ,:)], 'b - ');
title('H2')
axis equal

nexttile
plot ([ pflats3(1 ,:); pflate3(1 ,:)] , ...
[pflats3(2 ,:); pflate3(2 ,:)], 'b - ');
title('H3')
axis equal

nexttile
plot ([ pflats4(1 ,:); pflate4(1 ,:)] , ...
[pflats4(2 ,:); pflate4(2 ,:)], 'b - ');
title('H4')
axis equal

%% Computer Exercise 4
y = imread('compEx4im1.JPG');
y2 = imread('compEx4im2.JPG');
figure
% imagesc(y2)
hold on;

load('compEx4.mat');
P1 = K*[R1 t1];
P2 = K*[R2 t2];

% Viewing direction (principal axis) is the third row of P
vdP1 = P1(3,1:3);
vdP2 = P2(3,1:3);

vdP1Norm = vdP1 ./ sqrt(vdP1(1)^2 + vdP1(2)^2 + vdP1(3)^2)
vdP2Norm = vdP2 ./ sqrt(vdP2(1)^2 + vdP2(2)^2 + vdP2(3)^2)

% Finding the camera coordinates with this or -R'*t
C1 = pflat(null(P1));
C2 = pflat(null(P2));

% c1=-R1'*t1;
% c2=-R2'*t2;

Up = pflat(U);

% Extracting the coordinates
Cent1 = C1(1:3 ,1);
Cent2 = C2(1:3 ,1);

% plotting the U points, camera centers and the viewing direction
plot3(Up(1,:), Up(2,:), Up(3,:), '.', 'Markersize', 2);
hold on

plot3(Cent1(1,:), Cent1(2,:), Cent1(3,:), '.', 'color', 'red');
plot3(Cent2(1,:), Cent2(2,:), Cent2(3,:), '.', 'color', 'red');

quiver3(Cent1(1,:), Cent1(2,:), Cent1(3,:), vdP1(:,1), vdP1(:,2), vdP1(:,3), 10);
quiver3(Cent2(1,:), Cent2(2,:), Cent2(3,:), vdP2(:,1), vdP2(:,2), vdP2(:,3), 10);

%% plotting as 2d image
figure;
t = tiledlayout(1,2);
nexttile
hold on
up1 = pflat(P1*U);
colormap gray
imagesc(y)
plot(up1(1,:), up1(2,:), '.', 'Markersize', 2);

nexttile
hold on
up2 = pflat(P2*U);
colormap gray
imagesc(y)
plot(up2(1,:), up2(2,:), '.', 'Markersize', 2);
