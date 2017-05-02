function M2 = matrix_histogram2(samples,deltaq) 

    Bin_Width  = deltaq;
    num_measurements = size(samples, 1);
    num_angles = num_measurements/1000;
%     angles = pi*(0:num_angles-1)/num_angles;
%     angles = repmat(angles,(num_measurements/num_angles),(num_measurements/num_angles));
%     angles = angles(1:(num_measurements/num_angles)*num_angles)';
%     A = angles;

    H = zeros(num_measurements/num_angles, num_angles);
%     N = zeros(num_measurements/num_angles, num_angles);
    M = zeros(1,3);
    
    for i=1:num_angles;
        
        angle = samples(i,1); 
        H(:,i) = samples((i:num_angles:end),2); 
        [N,edges] = histcounts(H(:,i), 'BinWidth', Bin_Width);
        d = diff(edges)/2;
        centers = edges(1:end-1)+d;
        C = centers';  
        MA = [repmat(angle,length(N),1), C, N'];
        M = [M; MA];
        
        
    end
    M2 =M(2:end,:);
