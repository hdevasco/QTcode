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

method = {'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'};
num_measurements = size(samples, 1);
H = zeros(num_measurements/numAngles, numAngles);
M2 = zeros(1,3);
Bin_Width_Scott = zeros(num_measurements/(num_measurements/numAngles), 1);

if option > 8,
    
    num_bins = option;
    N = zeros(num_bins, numAngles);
    angles = pi*(0:numAngles-1)/numAngles;
    angles = repmat(angles,num_bins,(num_measurements/numAngles));
    angles = angles(1:num_bins*numAngles)';
    A = angles;
    
    for i=1:numAngles;
        
        H(:,i) = samples((i:numAngles:end),2);
        [N(:,i),edges] = histcounts(H(:,i), num_bins);
        d = diff(edges)/2;
        centers = edges(1:end-1)+d;
        C(:,i) = centers';
            A2 = A(:);
            C2 = C(:);
            N2 = N(:);
            M = [A2, C2, N2];
            
    end
    ind = find(M(:,3) > 0);
    M = M(ind,:);
end

for i=1:numAngles;
    angle = samples(i,1);
    H(:,i) = samples((i:numAngles:end),2);
    
    if (option ==8),
        
        Bin_Width_Scott(i) = 3.5*std(H(:,i))*((num_measurements/numAngles)^(-1/3));
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width_Scott(i));
        
    elseif (option ==7),
        
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width);
        
    elseif (option > 0) && (option < 7),
        
        [N,edges] = histcounts(H(:,i),'BinMethod', method{option});
        
        d = diff(edges)/2;
        centers = edges(1:end-1)+d;
        C = centers';
        MA = [repmat(angle,length(N),1), C, N'];
        M2 = [M2; MA];
    end
end
M =M2(2:end,:);
ind = find(M(:,3) > 0);
M = M(ind,:);

if strcmp(H_operator,'integral'),
    M_integral = zeros(1,4);
    for i=1:numAngles;
        angle = samples(i,1);
        N = M(:,3);
        centers = M(:,2);
        d = diff(centers)/2;
        edges = centers(1:end+1)-d;
        B = [edges(1:end-1)',edges(2:end)'];
        MA = [repmat(angle,length(N),1), B , N];
        M_integral = [M_integral; MA];
    end
    M =M_integral(2:end,:);
    ind = find(M(:,4) > 0);
    M = M(ind,:);
end
end
