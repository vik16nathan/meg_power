function output = convertTo3DArray(matrix, nrow, ncol)  
    % Initialize the 3D array
    output = zeros(nrow/3, ncol, 3);
    
    % Convert the matrix to as 3D array
    for n = 0:nrow/3-1
        output(n+1,:,1) = matrix(3*n+1,:,:);
        output(n+1,:,2) = matrix(3*n+2,:,:);
        output(n+1,:,3) = matrix(3*n+3,:,:);
    end
end
