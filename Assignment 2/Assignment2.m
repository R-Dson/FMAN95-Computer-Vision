% 2)
load('CompEx1data.mat');

% 
plot3(X(1,:), X(2,:), X(3,:), '.', 'Markersize', 2)
axis equal
hold on;
plotcams(P)

%% Computer Exercise 1.

i = 1;
y = imread(imfiles{i});

% X to camera i
px = P{i} * X;

% projection of px to image


% check if they are visible
visible = isfinite ( x{i}(1 ,:));

pf = pflat(px);

figure
axis equal
colormap gray
% Image
imagesc(y)
hold on;

% image points
plot( x { i }(1 , visible ) , x { i }(2 , visible ) , '*', 'Markersize', 3, 'color', 'red' );
% projection points
plot(pf(1, visible), pf(2, visible), '.', 'Markersize', 1, 'color', 'blue');

%% 
T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 1/10 1/10 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];

% applying H*X = H*T^-1*T*X
% new points are T * X
T1X = T1*X;
T2X = T2*X;

%Projecting to the image plane
T1Xp = pflat(T1X);
T2Xp = pflat(T2X);

% P * T^-1 gives a new camera. P{j}\T1 doesn't work?
for j = 1:size(P,2)
%     PT1{j} = P{j}\T1;
%     PT2{j} = P{j}\T2;
    PT1{j} = P{j}*inv(T1);
    PT2{j} = P{j}*inv(T2);
end

figure;
tiledlayout(1,2);
nexttile
plot3(T1Xp(1, :), T1Xp(2, :), T1Xp(3, :), '.');
hold on;
plotcams(PT1);
axis equal;

nexttile
plot3(T2Xp(1, :), T2Xp(2, :), T2Xp(3, :), '.');
hold on;
plotcams(PT2);
axis equal;

%%
y = imread(imfiles{i});

% X to camera i
px = PT1{i} * T1Xp;

% projection of px to image


% check if they are visible
visible = isfinite ( x{i}(1 ,:));

pf = pflat(px);

figure
axis equal
colormap gray

% Image
imagesc(y)
hold on;

% image points
plot( x { i }(1 , visible ) , x { i }(2 , visible ) , '*', 'Markersize', 3, 'color', 'red' );
% projection points
plot(pf(1, visible), pf(2, visible), '.', 'Markersize', 1, 'color', 'blue');


%% Computer Exercise 2.
format shortG
P1 = PT1{1};
[K1]=rq(P1);
K1 = K1./K1(3,3)

P2 = PT2{1};
[K2]=rq(P2);
K2 = K2./K2(3,3)

%% Computer Exercise 3
load('CompEx3data.mat');

% getting the x and y
xp1 = x{1}(1:2, :);
xp2 = x{2}(1:2, :);

% calculating the mean
xh1 = mean(xp1')';
xh2 = mean(xp2')';

% calculating the standard deviation
std1 = std(xp1')';
std2 = std(xp2')';

% N=[s 0 -sx^-; 0 s -sy^-; 0 0 1];, s is the standard deviation in our case
s1 = 1/std1(1);
s2 = 1/std1(2);
xhx = xh1(1);
xhy = xh1(2);
N1 = [s1 0 -s1*xhx; 0 s2 -s1*xhy; 0 0 1]

s1 = 1/std2(1);
s2 = 1/std2(2);
xhx = xh2(1);
xhy = xh2(2);
N2 = [s1 0 -s1*xhx; 0 s2 -s1*xhy; 0 0 1];
    
% using N*x to normalize x
N1x = N1*x{1};
N2x = N2*x{2};

figure;
tiledlayout(1,2)
nexttile
plot(N1x(1, :), N1x(2, :), '.', 'color', 'blue');
axis equal;

nexttile
plot(N2x(1, :), N2x(2, :), '.', 'color', 'blue');
axis equal;

%%
% using the normalized values N1x and N2x to build M

M1 = mBuild(Xmodel, N1x);
M2 = mBuild(Xmodel, N2x);

[u1,s1,v1] = svd(M1);
[u2,s2,v2] = svd(M2);
vsol1 = v1(:,end);
vsol2 = v2(:,end);

% using S^T*S to get the eigenvalues
e1 = s1'*s1;
e2 = s2'*s2;
smallest1 = e1(end, end);
smallest2 = e2(end, end);

M1v = norm(M1*vsol1);
M2v = norm(M2*vsol2);

% getting the vectors p_1, p_2 and p_3 to get camera matrix
P1 = [vsol1(1:4)'; vsol1(5:8)'; vsol1(9:12)'];
P2 = [vsol2(1:4)'; vsol2(5:8)'; vsol2(9:12)'];
% scale P back to before normalizing
P1 = N1\P1;
P2 = N2\P2;

P1X = P1*[Xmodel; ones(1, length(Xmodel))];
P2X = P2*[Xmodel; ones(1, length(Xmodel))];
p1x = pflat(P1X);
p2x = pflat(P2X);

im1 = imread("cube1.JPG");
im2 = imread("cube2.JPG");

figure;
tiledlayout(1,2);
nexttile
imagesc(im1);
% axis equal;
hold on;

plot(x{1}(1, :), x{1}(2, :), '+', 'color', 'red', 'MarkerSize', 10);
plot(p1x(1, :), p1x(2, :), '.', 'MarkerSize', 10, 'color', 'blue');

nexttile
imagesc(im2);
hold on;
% axis equal;
plot(x{2}(1, :), x{2}(2, :), '+', 'color', 'red', 'MarkerSize', 10);
plot(p2x(1, :), p2x(2, :), '.', 'MarkerSize', 10, 'color', 'blue');

% extracting model points
X = [Xmodel(1,startind); Xmodel(1,endind)];
Y = [Xmodel(2,startind); Xmodel(2,endind)];
Z = [Xmodel(3,startind); Xmodel(3,endind)];

% plotting the 3d model points and cameras P1 and P2
figure
plot3(X,Y,Z)
axis equal;
hold on;

% calculating camera center
C1 = -P1(1:3,1:3)\P1(:,4);
C2 = -P2(1:3,1:3)\P2(:,4);

plotcams({P1 P2});
plot3(C1(1), C1(2), C1(3), 'b.')
plot3(C2(1), C2(2), C2(3), 'b.')

% calculating rq
[r1, q1] = rq(P1);
[r2, q2] = rq(P2);
% dividing by r_i(3,3)
r1 = r1./r1(3,3)
r2 = r2./r2(3,3)
axis equal;

%% Computer Exercise 4.

im1 = imread("cube1.JPG");
im2 = imread("cube2.JPG");

[f1, d1] = vl_sift( single(rgb2gray(im1)), 'PeakThresh', 1);
[f2, d2] = vl_sift( single(rgb2gray(im2)), 'PeakThresh', 1);
figure
imagesc(im1);
hold on;
axis equal;
vl_plotframe(f1);

[matches ,scores] = vl_ubcmatch(d1, d2);

x1 = [f1(1,matches (1 ,:));f1(2,matches (1 ,:))];
x2 = [f2(1,matches (2 ,:));f2(2,matches (2 ,:))];

perm = randperm(size(matches ,2));
figure;
imagesc ([im1 im2]);
axis equal;
hold on;
plot([x1(1,perm (1:10)); x2(1,perm (1:10))+ size(im1 ,2)], ...
[x1(2,perm (1:10)); x2(2,perm (1:10))] ,'-');
hold off;

%% Computer Exercise 5.

for i = 1:size(x1,2)
    M = [[P1; P2] -[ [x1(:,i); 1] zeros(3,1); zeros(3,1) [x2(:,i); 1]]];
    [U,S,V] = svd(M);
    % getting the last column in V
    v = V(:, end);
    finalx(:, i) = v(1:4, 1);
end

% projection
xproj1 = pflat(P1*finalx);
xproj2 = pflat(P2*finalx);

% Plot the two images and points
figure;
tiledlayout(1,2);
nexttile
imagesc(im1);
hold on;
plot(x1(1,:), x1(2,:), '+', 'color', 'Y');
plot(xproj1(1,:), xproj1(2,:), '.', 'color', 'B');

nexttile
imagesc(im2);
hold on;
plot(x2(1,:), x2(2,:), '+', 'color', 'Y');
plot(xproj2(1,:), xproj2(2,:), '.', 'color', 'B');

% adding a threshhold of 3 pixels to the error points
good_points = (sqrt(sum((x1-xproj1(1:2,:)).^2)) < 3 & sqrt(sum((x2-xproj2(1:2,:)).^2)) < 3);
finalx = finalx(:, good_points );
xprojf = pflat(finalx);

% plot
figure;
plot3(xprojf(1, :), xprojf(2, :), xprojf(3, :), '.');
hold on;
plot3([Xmodel(1,startind); Xmodel(1,endind)], [Xmodel(2,startind); Xmodel(2,endind)], [Xmodel(3,startind); Xmodel(3,endind)], '-');
plotcams(P);