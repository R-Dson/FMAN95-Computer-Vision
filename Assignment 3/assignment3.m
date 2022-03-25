% 3)
clc
clear
load('compEx1data.mat');
img1 = imread('kronan1.JPG');
img2 = imread('kronan2.JPG');

% Computer Exercise 1.

N1 = normalizeM(x{1});
N2 = normalizeM(x{2});

% N1 = eye(3);
% N2 = eye(3);

% using N*x to normalize x
x1n = N1*x{1};
x2n = N2*x{2};

%setting up M for the 8 points
for i = 1:length(x{1})
    xx = x2n(:,i)*x1n(:,i)';
    M(i,:) = xx(:)';
end

% SVD solution
[U, S, V] = svd(M);
v = V(:, end); 

% both are small
sdiag = diag(S);
minst = min(sdiag)
absmv = norm(M*v)

% gets F tilde
Fn = reshape (v, [3 3]);
% solving min det(F)=0 ||F tidle - F||
[U, S, V] = svd(Fn);
% setting smallest signular value to 0 so we get F
S(end,end)=0;
Fn=U*S*V';
% now det(Fn) is -3.2610e-19
determinant = det(Fn)
xTFx = x2n(:,1)'*Fn*x1n(:,1)
% un-normalized
F = N2'*Fn*N1;
F = F./F(3, 3);
determinantF = det(F)

% computing the epipolar lines
l=F*x{1};
l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1]));

ran = randi([1 length(x{1})],1, 20);
points = x{2}(:, ran);
lpoints = l(:, ran);

% Plot
figure;
imagesc(img2);
hold on;

plot(points(1,:), points(2,:), '.', 'color' , 'g', 'MarkerSize',15);
rital(lpoints);
axis equal;

% corrected mine with the one you provided, both give the same
d = abs(sum(l.*x{2}));
meandistNew = mean(d)
meandist = sum(abs(sum(l .* x {2})))/length(mean(l .* x {2}))

figure;
plot ( diag ( x2n'* Fn * x1n ));

figure;
hist ( abs ( sum ( l .* x {2})) ,100);
%% Computer Exercise 2.
 
 
%P1 is always [I 0]
P1 = [eye(3) zeros(3, 1)];
 
% e2 (t part) is in the nullspace of F^T
e2 = null(F');
% A = e_2x * F
A = [0 -e2(3) e2(2); 
    e2(3) 0 -e2(1); 
    -e2(2) e2(1) 0] * F;
% P_2 = [A|e2]
P2 = [A e2];
 
% DLT algorithm
% normalize P1 and P2
NormP1 = N1*P1;
NormP2 = N2*P2;
% normalize the points 
nx1 = N1 * x{1};
nx2 = N2 * x{2};
 
for i = 1:length(x{1})
    msvd = [NormP1 -nx1(:, i) zeros(3, 1); 
            NormP2 zeros(3, 1) -nx2(:, i)];
    [U, S, V] = svd(msvd);
    xpoints(:, i) = V(1:4, end);
end
 
% switching columns 3 and 4 in camera and row 3 and 4 in 3d points
xpoints = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]*xpoints;
P1 = P1*[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
 
 
% using the pflat to project to the image plane
xp = pflat(P1*xpoints);
xpoints = pflat(xpoints);
 
% Plot
figure;
imagesc(img1);
hold on;
% image points
plot(x{1}(1,:), x{1}(2,:), '+', 'color', 'y', 'MarkerSize', 5);
hold on;
% traiangulated points
plot(xp(1,:), xp(2,:), '.', 'color', 'b', 'MarkerSize', 5);
axis equal;
 
% 3D points plot
figure;
plot3(xpoints(1, :), xpoints(2, :), xpoints(3, :), '.');
% axis equal

%% Computer Exercise 3.
clear
load('compEx1data.mat');
load('compEx3data.mat');

img1 = imread('kronan1.JPG');
img2 = imread('kronan2.JPG');

nx1 = K\x{1};
nx2 = K\x{2};

% implementing algorithm again
% no need to normalize now
for i = 1:length(x{1})
    xx = nx2(:,i)*nx1(:,i)';
    M(i,:) = xx(:)';
end

% solving min ||v||=1 for ||Mv||^2
[U, S, V] = svd(M);
v = V(:, end);
res = reshape(v, [3 3]);

[U, S, V] = svd(res);

% finding norm and minimum signular value
sdiag = diag(S);
minst = min(sdiag)
absmv = norm(M*v)

% making sure E has the singular values 1,1,0
if det(U*V') >0
    E = U * diag ([1 1 0])* V';
else
    V = -V;
    E = U * diag ([1 1 0])* (-V)';
end

E = E./E(3, 3)

% checking some values
xTEx = nx2(:,1)' * E * nx1(:,1)
xTEx = nx2(:,2)' * E * nx1(:,2)
xTEx = nx2(:,3)' * E * nx1(:,3)

% getting the F matrix using  K'^-1 E K^-1 = F
F = K'\E/K;

% computing the epipolar lines like before
l = F*x{1};
l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1]));

ran = randi([1 length(x{1})],1, 20);
pointsnew = x{2}(:, ran);
lpointsnew = l(:, ran);

% Plot
figure;
imagesc(img1);
hold on;

% image points
plot(pointsnew(1,:), pointsnew(2,:), '.', 'color', 'green', 'MarkerSize', 12);
rital(lpointsnew);
axis equal;

figure;
hist ( abs ( sum ( l .* x {2})) ,100);
%% Computer Exercise 4.
% U = [1/sqrt(2) -1/sqrt(2) 0; 1/sqrt(2) 1/sqrt(2) 0; 0 0 1];
% V = [1 0 0; 0 0 -1; 0 1 0];
W = [0 -1 0; 1 0 0; 0 0 1];
u3 = U(:,3);

% using U and V from before to get P
P0=[eye(3) zeros(3,1)];
P1=[U*W*V' u3];
P2=[U*W*V' -u3];
P3=[U*W'*V' u3];
P4=[U*W'*V' -u3];
Ps = {P1 P2 P3 P4};
% Updated code
% Using my code from assignment 4
bestNumP = 0;
bestP=[];
bestxpoints3d=[];
bestxpoints2d=[];
numbestxpoints = 0;
ind = 0;
for i = 1:length(Ps)
    
    xpoints = svdCalc(P0, Ps{i}, nx1, nx2);
    
    P0xpoints = P0*xpoints;
    Pixpoints = Ps{i}*xpoints;
    
    points = sum(Pixpoints(3,:)>0) + sum(P0xpoints(3,:)>0);
    
    if points > bestNumP
        bestxpoints3d = xpoints;
        bestxpoints2d = Pixpoints;
        bestP = Ps{i};
        bestNumP = points;
        ind = i;
    end
end


% Plot the best solution
figure;
plot3(bestxpoints3d(1,:), bestxpoints3d(2,:), bestxpoints3d(3,:), '.', 'MarkerSize', 5);
hold on;
plotcams({P0, bestP});
axis equal;

% undo normalize
best = pflat(K*bestxpoints2d);

% Plot the image, points and projections
figure;
imagesc(img2);
hold on;
plot(x{2}(1,:), x{2}(2,:), '+', 'color', 'y', 'MarkerSize', 5);
plot(best(1,:), best(2,:), '.', 'color', 'b', 'MarkerSize', 5);
axis equal;


%%
% W = [0 -1 0; 1 0 0; 0 0 1];
% u3 = U(:,3);
% 
% % using U and V from before to get P
% P0=[eye(3) zeros(3,1)];
% P1=[U*W*V' u3];
% P2=[U*W*V' -u3];
% P3=[U*W'*V' u3];
% P4=[U*W'*V' -u3];
% Ps = {P1 P2 P3 P4};
% % In hindisght, this could probably be done in a forloop
% % doing the svd to find the points from V
% xpoints1 = svdCalc(P0, P1, nx1, nx2);
% xpoints2 = svdCalc(P0, P2, nx1, nx2);
% xpoints3 = svdCalc(P0, P3, nx1, nx2);
% xpoints4 = svdCalc(P0, P4, nx1, nx2);
% 
% % need to pflat to get the 3d projected points
% 
% P0xpoints1 = P0*xpoints1;
% P0xpoints2 = P0*xpoints2;
% P0xpoints3 = P0*xpoints3;
% P0xpoints4 = P0*xpoints4;
% 
% P1xpoints1 = P1*xpoints1;
% P2xpoints2 = P2*xpoints2;
% P3xpoints3 = P3*xpoints3;
% P4xpoints4 = P4*xpoints4;
% 
% % z needs to be in front of the camera (>0)
% 
% P1xp1pos = sum(P1xpoints1(3,:)>0) + sum(P0xpoints1(3,:)>0);
% P2xp2pos = sum(P2xpoints2(3,:)>0) + sum(P0xpoints2(3,:)>0);
% P3xp3pos = sum(P3xpoints3(3,:)>0) + sum(P0xpoints3(3,:)>0);
% P4xp4pos = sum(P4xpoints4(3,:)>0) + sum(P0xpoints4(3,:)>0);
% 
% %finding the best camera matrix
% [M,I] = max([P1xp1pos, P2xp2pos, P3xp3pos, P4xp4pos])
% % getting it
% best = {P1xpoints1, P2xpoints2, P3xpoints3, P4xpoints4};
% best3d = {xpoints1, xpoints2, xpoints3, xpoints4};
% best = best{I};
% best3d = pflat(best3d{I});
% 
% % Plot the best solution
% figure;
% plot3(best3d(1,:), best3d(2,:), best3d(3,:), '.', 'MarkerSize', 5);
% hold on;
% Ps={P1,P2,P3,P4};
% plotcams({P0, Ps{I}});
% axis equal;
% 
% % undo normalize
% best = pflat(K*best);
% 
% % Plot the image, points and projections
% figure;
% imagesc(img2);
% hold on;
% plot(x{2}(1,:), x{2}(2,:), '+', 'color', 'y', 'MarkerSize', 5);
% plot(best(1,:), best(2,:), '.', 'color', 'b', 'MarkerSize', 5);
% axis equal;

function xpoints = svdCalc(P0, P, x1, x2)
    for i = 1:length(x1)
        msvd = [P0 -x1(:, i) zeros(3, 1); 
                P zeros(3, 1) -x2(:, i)];
        [U, S, V] = svd(msvd);
        v = V(1:4, end);
        xpoints(:, i) = -pflat(v);
    end
end

function N = normalizeM(M)
    tM = M';
    mM = mean(tM)';
    stdM = std(tM)';
    % N=[s 0 -sx^-; 0 s -sy^-; 0 0 1];, s is the standard deviation in our case
    s1 = 1/stdM(1);
    s2 = 1/stdM(2);
    N = [s1 0 -s1*mM(1); 
        0 s2 -s2*mM(2); 
        0 0 1];
end
