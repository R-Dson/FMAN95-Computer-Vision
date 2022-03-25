%% Computer Exercise 1.
clc
clear
ima = imread('a.jpg');
imb = imread('b.jpg');

imagesc(ima);
imagesc(imb);
axis equal;

% Sift features
[ fA dA ] = vl_sift ( single ( rgb2gray ( ima )) );
[ fB dB ] = vl_sift ( single ( rgb2gray ( imb )) );

% Matches
[matches ,scores] = vl_ubcmatch(dA, dB);

% Getting matches
xA = [fA(1:2 , matches (1 ,:)); ones(1, length(matches)) ];
xB = [fB(1:2 , matches (2 ,:)); ones(1, length(matches)) ];

bestInlier = 0;
bestH = [];
minP = 4; % min points
maxD = 5; % threshold

% RANSAC algorithm using 4 points, 5 distance threshhold and 
% log(0.02)/log(1-0.9^4) rounded up is 4 as N. Using 10 instead 
% for good measure.
for i = 1:10
    % selecting random points
    rpoints=zeros(1,minP);
    for j = 1:minP
        rpoints(j) = randi([1 length(xA)]); 
    end
    rxA = xA(:, rpoints);
    rxB = xB(:, rpoints);
    
    % Calculating M, SVD and getting last row of V
    M = mBuild(rxA, rxB);
    [U, S, V] = svd(M);
    v = V(:, end);
    H = reshape(v(1:9), [3,3])';
    
    % homography projection
    pnts = H*xA;
    
    numInliers = 0;
    for j = 1:size(xA, 2)
        pt1 = pnts(:,j);
        pt1 = pt1/pt1(3);
        pt2 = xB(:, j);
        
        d = norm(pt1 - pt2);
        if d < maxD
            numInliers = numInliers + 1;
        end
    end
    % finds the most inliers
    if  bestInlier < numInliers
        bestInlier = numInliers;
        bestH = H;
    end
end

% given in the assignment
Htform = projective2d ( bestH');
% Creates a transfomation that matlab can use for images .
% Note : MATLAB uses transposed transformations .

Rout = imref2d ( size ( ima ) ,[ -200 800] ,[ -400 600]);
% Sets the size and output bounds of the new image .

[ Atransf ] = imwarp (ima , Htform , 'OutputView' , Rout );
% TransfoRMS the image

Idtform = projective2d ( eye (3));
[ Btransf ] = imwarp (imb , Idtform , 'OutputView' , Rout );
% Creates a larger version of the second image

AB = Btransf ;
AB ( Btransf < Atransf ) = Atransf ( Btransf < Atransf );
% Writes both images in the new image . %( A somewhat hacky solution is needed
% since pixels outside the valid image area are not always zero ...)

imagesc ( Rout . XWorldLimits , Rout . YWorldLimits , AB );
% Plots the new image with the correct axes

axis equal


%% Computer Exercise 2.
load('compEx2data.mat');
im1 = imread('im1.jpg');
im2 = imread('im2.jpg');

x1 = [x{1}; ones(1, length(x{1}) )];
x2 = [x{2}; ones(1, length(x{2}) )];
%normalize
x1n = K\x1;
x2n = K\x2;

% variables we need
bestInlier = 0;
bestInliers = [];
bestE = [];
minP = 5; % min points
maxD = 5; % threshold
% RANSAC to find best estimation
for i = 1:10
    % selecting random points
    rpoints=zeros(1,minP);
    for j = 1:minP
        rpoints(j) = randi([1 length(x1)]); 
    end
    
    % 5 points solver for random points to get list of E's
    Elist = fivepoint_solver(x1n(:, rpoints), x2n(:, rpoints));
    
    for j = 1:size(Elist, 2)
        % temp variables
        numInliers = 0;
        inliers = [];
        % similar to assignment 3
        
        % solving for F
        F = (K')\Elist{j}/K;
        
        for k = 1:length(x1)
            % sill similar to assignment 3
            l=F*x1(:,k);
            l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1]));
            % distance to epipolarline
            d = abs ( sum ( l .* x2(:,k)));
            
            if d < maxD
                numInliers = numInliers + 1;
                inliers(end+1) = k;
            end
        end
        % finds the most inliers
        if  bestInlier < numInliers
            bestInlier = numInliers;
            bestE = Elist{j};
            bestInliers = inliers;
        end
    end
end

% getting inliers. this format is what ComputeReprojectionError uses
u = {x1(:, bestInliers), x2(:, bestInliers)};


% calculating SVD
% [U, S, V] = svd(bestE./bestE(3,3));

% UPDATE: Also added this if check to determine if V should be positive or
% negative.
[U, S, V] = svd(bestE./bestE(3,3));
W = [0 -1 0; 1 0 0; 0 0 1];
u3 = U(:,3);

if det(U * V') < 0
    V = -V;
end

% doing like we did in assignment 3
% unnormalized camera
P0=K*[eye(3) zeros(3,1)];
P1=K*[U*W*V' u3];
P2=K*[U*W*V' -u3];
P3=K*[U*W'*V' u3];
P4=K*[U*W'*V' -u3];
Ps = {P1 P2 P3 P4};

bestNumP = 0;
bestP=[];
bestxpoints=[];
for i = 1:length(Ps)
%     points = 0;
    
    xpoints = svdCalc(P0, Ps{i}, u{1}, u{2});
    
    % I tried to use this at first similar to assignment 3
    % but it's not always correct and sometimes gives a
    % flat reconstruction. 
    % Update: This has now been fixed by checking if det(U * V') < 0
    
    % Update: Using this again after comment and rewrote it.
    % calculating number of points in front
    P0xpoints = P0*xpoints;
    Pixpoints = Ps{i}*xpoints;
    p0front = (Pixpoints(3,:)>0);
    p1front = (P0xpoints(3,:)>0);
    points = sum( (p0front + p1front ) == 2);
    
    % Using depth gives the correct reconstruction
    % calculating depth
%     d1 = depthCalc(P0, xpoints);
%     d2 = depthCalc(Ps{i}, xpoints);
%     s1 = sum(d1 > 0);
%     s2 = sum(d2 > 0);

%     if s1 > 0 && s2 > 0
%         points = points + 1;
%     end
    if points > bestNumP
        bestxpoints = xpoints;
        bestP = Ps{i};
        bestNumP = points;
    end
end

% section 4 contains this code that I used and they are given in the assignment
P = {P0, bestP};
[err, res] = ComputeReprojectionError(P, bestxpoints, u);
figure;
hist ( res, 100 );
% RMS = sqrt((sum of n_i ) / number of n)
RMS = sqrt(sum(res.^2)/length(res))

% The reconstruction
figure;
plot3(bestxpoints(1, :), bestxpoints(2, :), bestxpoints(3, :), '.', 'color', 'b');
hold on;
% plotcams(P);
axis equal;

%% Computer Exercise 3.
% using data from Comp. Ex. 2
% some value for gamma to start at
gammak = 0.5;
% err2 will get overwritten later so start values doesnt matter
err2 = 0;
% bestxpoints is our U
U = bestxpoints;

errors = [err];
for i = 1:10
    % init values. parameters (P, U, u)
    [r, J] = LinearizeReprojErr(P, U, u);
    
    % calculating deltav
    deltav = - gammak *J'* r ;
    % parameters (deltav, P, U)
    [ Pnew , Unew ] = update_solution ( deltav ,P , U );
    % parameters (P, U, u)
    [err2, res] = ComputeReprojectionError(Pnew, Unew, u);
    % we now have a error difference of err and err2
    errors(i+1) = err2;
    %repeat until we get a error less than previous for minimizing error
    while err < err2
        gammak = gammak/3;
        % some other gammak I tried
%         gammak = gammak/5;
%         gammak = gammak/10;
%         gammak = gammak*99/100; % this one is interesting
%         gammak = gammak*exp(-1.5);

        deltav = - gammak *J'* r ;
        % parameters (deltav, P, U)
        [ Pnew , Unew ] = update_solution ( deltav ,P , U );
        % parameters (P, U, u)
        [err2, res] = ComputeReprojectionError(Pnew, Unew, u);
        
        errors(i+1) = err2;
    end
    err = err2;
    % saving the new values with smaller error
    P = Pnew;
    U = Unew;
%     bestxpoints = Unew;
end
figure;
% including the one from comp. ex. 2 too so we start at 0 with 11 points
plot(linspace(0,10,11), errors);

RMS = sqrt(sum(res.^2)/length(res))



%% Computer Exercise 4.
% using data from Comp. Ex. 2 (run as section) and skip comp ex 3

% δv = −(J(vk)T J(vk) + λI)−1J(vk)T r(vk)
lambda = 1/10;
err2 = 0;
% similar to comp ex 3
errors = [err];
U = bestxpoints;

for i = 1:10
    % Computes the LM update .
    [r, J] = LinearizeReprojErr(P, U, u);
    % the (J(vk)T J(vk) + λI)−1 part
    C = J'* J + lambda * speye ( size (J ,2));
    % the J(vk)T r(vk) part
    c = J'* r ;
    % we get δv = −(J(vk)T J(vk) + λI)−1J(vk)T r(vk)
    deltav = -C \ c ;

    [ Pnew , Unew ] = update_solution ( deltav ,P , U );
    [err2, res] = ComputeReprojectionError(Pnew, Unew, u);
    P = Pnew;
    errors(i+1) = err2;
    U= Unew;
end
plot(linspace(0,10,11), errors);
title('lambda = 1/10')
RMS = sqrt(sum(res.^2)/length(res))






% function xdepth = depthCalc(P, x)
%     % depth = sign(deth(A))/||A|| * [A_3' a_3]X
%     A = P(:,1:3);
%     A3 = A(3,:);
%     a = P(:,4);
%     a3 = a(3,1);
%     % not sure why A3 works and not A3'
%     xdepth = (sign(det(A))/norm(A3)) * [A3 a3]*x;
% end

function xpoints = svdCalc(P0, P, x1, x2)
    for i = 1:length(x1)
        msvd = [P0 -x1(:, i) zeros(3, 1); 
                P zeros(3, 1) -x2(:, i)];
        [U, S, V] = svd(msvd);
        v = V(:, end);
        xpoints(:, i) = pflat(v(1:4, 1));
    end
end
