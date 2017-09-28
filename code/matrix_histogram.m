function M = matrix_histogram(samples, option, deltaq) 
% M = MATRIX_HISTOGRAM(samples, option) 
%  This functions makes a matrix histogram using as input samples a matrix
%  Mx3 of columns [angles, measurements of quadratures, measures number in each bin]
%  and a option using this pattern:
%   option = 
%       1 = 'auto'
%       2 = 'scott'
%       3 = 'fd'
%       4 = 'integers'
%       5 = 'sturges'
%       6 = 'sqrt'
%       7 = 'BinWidth'(Option 7 specifies the width of the bin, which must be initialized 
%       in matrix_histogram(sample, 7, bin width(scalar))).
%       values grater than 7 indicates a number of bins.

num_measurements = size(samples, 1);
num_angles = num_measurements/1000;
method = {'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'};

% Specify number of bins
if option > 7,
    num_bins = option;
    angles = pi*(0:num_angles-1)/num_angles;
    angles = repmat(angles,num_bins,(num_measurements/num_angles));
    angles = angles(1:num_bins*num_angles)';
    A = angles;

    H = zeros(num_measurements/num_angles, num_angles);
    N = zeros(num_bins, num_angles);

    for i=1:num_angles;
         
        H(:,i) = samples((i:num_angles:end),2); 
        [N(:,i),edges] = histcounts(H(:,i), num_bins);
        d = diff(edges)/2;
        centers = edges(1:end-1)+d;
        C(:,i) = centers';        
    end
    
    A2 = A(:);
    C2 = C(:);
    N2 = N(:);
    
    ind = find(N2 > 0);
    
    A3 = A2(ind);
    C3 = C2(ind);
    N3 = N2(ind);
    
    M = [A3, C3, N3];
else
    % Specified BinWidth
    if (option ==7),
    Bin_Width  = deltaq;
    num_measurements = size(samples, 1);
    num_angles = num_measurements/1000;
    
    H = zeros(num_measurements/num_angles, num_angles);
    M = zeros(1,4);
    % The matrix B is the matrix N by 2 that represents in each line the limits of the bin.
    B = zeros(length(edges-1),2);
    for i=1:num_angles;
        
        angle = samples(i,1); 
        H(:,i) = samples((i:num_angles:end),2); 
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width);
        B = [edges(1:end-1)',edges(2:end)'];
        MA = [repmat([angle,length(N),1), B , N'];
%         d = diff(edges)/2;
%         centers = edges(1:end-1)+d;
%         C = centers';  
%         MA = [repmat(angle,length(N),1), C, N'];
         M = [M; MA];
        
        end
    end
    M =M(2:end,:);
else
    % Specify method
    if (option > 0) && (option < 7),

        H = zeros(num_measurements/num_angles, num_angles);
        M = zeros(1, 3);

        for i=1:num_angles;

                angle = samples(i,1);
                H(:,i) = samples((i:num_angles:end),2); 
                [N,edges] = histcounts(H(:,i),'BinMethod', method{option});
                d = diff(edges)/2;
                centers = edges(1:end-1)+d;
                C = centers';
                MA = [repmat(angle,length(N),1), C, N'];
                M = [M; MA];

         end

             M =M(2:end,:);
    
        end
    end
end 

end