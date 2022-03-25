function M = mBuild(X, coord)
    % height: 3*number of coords. 
    % width: 12 from X_i (4x1) and two zeros vectors (4x1), so (4x3)
    num = size(coord, 2);
    num1 = size(coord, 1);
    M = zeros(3 * num, 3*num1 + num);
    X = X';
    
    % building the M matrix [x_1 0 0 -x_1 0 0...; 0 x_1 0 -y_1 0 0...; 0 0
    % x_1 -1 0 0 0...] and so on
    j = 1;
    for i = 1:num
        M(j, 1:3) = X(i,:);
        M(j+1, 4:6) = X(i,:);
        M(j+2, 7:9) = X(i,:);
        
        M(j, 3*num1+i) = -coord(1,i);
        M(j+1, 3*num1+i) = -coord(2,i);
        M(j+2, 3*num1+i) = -coord(3,i);
        
        j = j + 3;
    end
end

