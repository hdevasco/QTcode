function M = matrix_histogram_new(numAngles, samples, option, H_operator,deltaq)
%  matrix_histogram.m returns the matrix of measurements through the use of the histogram in the data
%  This functions makes a matrix histogram using as input samples a matrix
%  Mx3 of columns [angles, measurements of quadratures,number of
%  observations] and return an array with [angle, quadrature value of the center of the box, number of
%  counts am each bin] for the input H_operator = center and the array with [angle, left edge of bin, right
%  edge of bin, number of counts am each bin] for the H_operator = integral input.
%  the "H_operator = center" option will allow homodyne_loss_measurement.m to calculate the measurement
%  operator in the center of the bin while "H_operator = integral" will allow coarse_measurement.m to
%  calculate the integrated measurement operator along the bin.
%  The option option allows you to choose how the histogram is constructed either by setting the width
%  of the bin directly or by a method of optimal width of the ones present in the statistic.
%   option =
%       1 = 'auto'
%       2 = 'scott'
%       3 = 'fd'
%       4 = 'integers'
%       5 = 'sturges'
%       6 = 'sqrt'
%       7 = 'BinWidth'(Option 7 specifies the width of the bin, which must be initialized
%       in matrix_histogram(sample, 7, bin width(scalar))).
%       8 = Determines the actual width of the Scott method without rounding to two significant
%       decimal places for the case of integration into the edges of the box.
%       values grater than 8 indicates a number of bins.

num_measurements = size(samples, 1);
method = {'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'};

% Specify number of bins
if option > 8,
    
    num_bins = option;
    angles = pi*(0:numAngles-1)/numAngles;
    angles = repmat(angles,num_bins,(num_measurements/numAngles));
    angles = angles(1:num_bins*numAngles)';
    
    A = angles;
    H = zeros(num_measurements/numAngles, numAngles);
    N = zeros(num_bins, numAngles);
    M2 = zeros(1,4);
    
    for i=1:numAngles;
        angle  = samples(i,1);
        H(:,i) = samples((i:numAngles:end),2);
        [N(:,i),edges] = histcounts(H(:,i), num_bins);
        
        if strcmp(H_operator,'center')
            
            d = diff(edges)/2;
            centers = edges(1:end-1)+d;
            C(:,i) = centers';
            
            A2 = A(:);
            C2 = C(:);
            N2 = N(:);
            
            ind = find(N2 > 0);
            
            A3 = A2(ind);
            C3 = C2(ind);
            N3 = N2(ind);
            
            M = [A3, C3, N3];
            
        elseif strcmp(H_operator,'integral')
            
            B = zeros(length(edges-1),2);
            B = [edges(1:end-1)',edges(2:end)'];
            MA = [repmat(angle,length(N(:,i)),1), B , N(:,i)];
            M2 = [M2; MA];
            
            M =M2(2:end,:);
        end
    end
end

if (option ==8),
    
    H = zeros(num_measurements/numAngles, numAngles);
    M = zeros(1,3);
    M2 = zeros(1,4);
    Bin_Width_Scott = zeros(num_measurements/(num_measurements/numAngles), 1);
    for i=1:numAngles;
        angle1 = samples(i,1);
        H(:,i) = samples((i:numAngles:end),2);
        Bin_Width_Scott(i) = 3.5*std(H(:,i))*((num_measurements/numAngles)^(-1/3));
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width_Scott(i));
        
        if strcmp(H_operator,'center'),
            
            d = diff(edges)/2;
            centers = edges(1:end-1)+d;
            C = centers';
            MA = [repmat(angle1,length(N),1), C, N'];
            M = [M; MA];
            M =M(2:end,:);
            
        elseif strcmp(H_operator,'integral'),
            
            B = zeros(length(edges-1),2);
            B = [edges(1:end-1)',edges(2:end)'];
            MA = [repmat(angle1,length(N),1), B , N'];
            M2 = [M2; MA];
            M =M2(2:end,:);
        end
    end
end

if (option ==7),
    
    Bin_Width  = deltaq;
    H = zeros(num_measurements/numAngles, numAngles);
    M = zeros(1,3);
    M2 = zeros(1,4);
    
    for i=1:numAngles;
        angle2 = samples(i,1);
        H(:,i) = samples((i:numAngles:end),2);
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width);
        
        if strcmp(H_operator,'center'),
            
            d = diff(edges)/2;
            centers = edges(1:end-1)+d;
            C = centers';
            MA = [repmat(angle2,length(N),1), C, N'];
            M = [M; MA];
            M =M(2:end,:);
            
        elseif strcmp(H_operator,'integral'),
            
            B = zeros(length(edges-1),2);
            B = [edges(1:end-1)',edges(2:end)'];
            MA = [repmat(angle2,length(N),1), B , N'];
            M2 = [M2; MA];
            M =M2(2:end,:);
            
        end
    end
end

if (option > 0) && (option < 7),
    
    H = zeros(num_measurements/numAngles, numAngles);
    M = zeros(1, 3);
    M2 = zeros(1,4);
    
    for i=1:numAngles;
        
        angle3 = samples(i,1);
        H(:,i) = samples((i:numAngles:end),2);
        [N,edges] = histcounts(H(:,i),'BinMethod', method{option});
        
        if strcmp(H_operator,'center'),
            
            d = diff(edges)/2;
            centers = edges(1:end-1)+d;
            C = centers';
            MA = [repmat(angle3,length(N),1), C, N'];
            M = [M; MA];
            M =M(2:end,:);
            
        elseif strcmp(H_operator,'integral'),
            
            B = zeros(length(edges-1),2);
            B = [edges(1:end-1)',edges(2:end)'];
            MA = [repmat(angle3,length(N),1), B , N'];
            M2 = [M2; MA];
            M =M2(2:end,:);
            
        end
    end
end


