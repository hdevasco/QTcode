function M = matrix_histogram_new(numAngles, samples, option, H_operator,deltaq)

% matrix_histogram_new constructs the matrix of histograms M according to the sample and the measurement operator that you 
% want to use for the reconstruction of the state: measurement operator representing the measurement in 
% the center of the box (center) and measuring operator along the length of the box ( integral). 
% If the option chosen for the use of the measuring operator is "center", M will be an array with columns 
% (angle, measured in the center of the box, number of counts of the box). If the option chosen for the use 
% of the measurement operator is "integral", M will be a matrix with columns (angle, left edge of the box, 
% right edge of the box, number of counts of the box).

% Constructs the histogram matrix using the edges of the boxes whose width is chosen by the Leonhardt method 
% using the measurement operator along the box.

if nargin == 2
  option = 'bin_leonhardt';
  H_operator = 'integral';
end

num_measurements = size(samples, 1);

angles = samples(1:numAngles);

% constructs the histograms using the homodyne measurements of the sample by choosing the method to 
% calculate the width of the box.

for i=1:numAngles
    QuadHist(i).angle = angles(i);
    QuadHist(i).allQuads = samples((i:numAngles:end),2);
    if strcmp(option,'number_bins')
        num_of_bins = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, num_of_bins);
    elseif strcmp(option,'scott_true')
        Bin_Width_Scott = 3.5*std(QuadHist(i).allQuads)*((num_measurements/numAngles)^(-1/3));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width_Scott);
    elseif strcmp(option,'bin_width')
        Bin_Width = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
    elseif strcmp(option,'bin_leonhardt')
        n = n_quadrature(samples,num_measurements);
        Bin_Width = pi/(2*sqrt(2*n+1));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
    elseif any(strcmp(option,{'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'}))
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinMethod', option);
     else 
         error('Error when the option is not (number_bins, bin_width,bin_leonhardt, scott_true, auto, scott, fd, integers, sturges or sqrt)')
    end
end

% Constructs M so that the measurement operator represents the measurement that occurs exactly in the center 
% of the histogram box.

if strcmp(H_operator,'center')
    
    for i = 1:numAngles
        d = diff(QuadHist(i).edges)/2;
        QuadHist(i).centers = QuadHist(i).edges(1:end-1)+d;
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1), QuadHist(i).centers.',QuadHist(i).counts.'];
    end
   
% Constructs M so that the measurement operator represents the measure that occurs along the length of the 
% histogram box.    

elseif strcmp(H_operator,'integral')
    
    for i = 1:numAngles
        QuadHist(i).edgesArray = zeros(length(QuadHist(i).edges-1),2);
        QuadHist(i).edgesArray= [QuadHist(i).edges(1:end-1).',QuadHist(i).edges(2:end).'];
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1) , QuadHist(i).edgesArray ,QuadHist(i).counts.'];
    end
  
end

M = vertcat(QuadHist.M);
M = M(M(:,end)> 0,:);

end
