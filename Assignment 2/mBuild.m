function M = mBuild(X, coord)
    % height: 3*number of coords. 
    % width: 12 from X_i (4x1) and two zeros vectors (4x1), so (4x3)
    num = size(coord, 2);
    M = zeros(3 * num, 12 + num);
    X = [X' ones(num, 1)];
    
    % building the M matrix [x_1 0 0 -x_1 0 0...; 0 x_1 0 -y_1 0 0...; 0 0
    % x_1 -1 0 0 0...] and so on
    j = 1;
    for i = 1:num
        M(j, 1:4) = X(i,:);
        M(j+1, 5:8) = X(i,:);
        M(j+2, 9:12) = X(i,:);
        
        M(j, 12+i) = coord(1,i);
        M(j+1, 12+i) = coord(2,i);
        M(j+2, 12+i) = coord(3,i);
        
        j = j + 3;
    end
end

