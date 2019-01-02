function M = matrix_histogram(numAngles, samples, option, H_operator, ...
                              deltaq)

% matrix_histogram discretizes the continuous values of homodyne measurements
% giving as result a matrix containing the discretized values.
%   Required inputs:
%      numAngles = number of evenly spaced phases from 0 to Pi. 
%                  In https://arxiv.org/pdf/1805.07414.pdf this is 'm' (= 20).
%      samples = matrix containing the results of the homodyne
%                measurements, which can be produced by the function homodyne_samples.
%   Optional inputs:
%      option = indicates the method used to set the width of the histogram bin:
%            -If 'number_of_bins' is chosen, the histogram will be constructed using 
%             a number deltaq of bins.
%            -If 'bin_width' is chosen, the histogram will be constructed using bins 
%             with width given by deltaq. 
%            -If 'bin_leonhardt' is chosen, the bin width will be
%             given by Leohnardt's equation, published in [Leonhardt,
%             Munroe, Kiss, Richter, Raymer "Sampling of photon
%             statistics and density matrix using homodyne
%             detection" Optics Communications 127, 144 (1996)],
%             Eq. (41)
%            -If 'scott_true' is chosen, the width will be given by
%             Scott's formula from [Scott "Scott’s rule" WIREs
%             Comp. Stat. 2, 497 (2010)], Eq (1).  Matlab provides
%             an implementation of Scott's formula, but Matlab's
%             implementation applies some rounding.
%            -The other options ('auto', 'scott', 'fd', 'inteiros', 'sturges', 'sqrt') use 
%             the automatic binning algorithm from the Matlab function histcounts. This 
%             binning algorithm returns bins with a uniform width, chosen to cover the 
%             range of elements in the distribution and reveal the underlying shape of the
%             distribution.
%      H_operator = 'center' chooses 3 column output
%                    format. 'integral' chooses 4 column output format. The 3 and
%                    4 column outputs are interpreted differently by
%                    make_measurement_struct.m, which will compute POVM elements
%                    as if the measurements occurred at the center of each
%                    histogram bin when given 3 columns, and will integrate the
%                    POVM elements over each bin when given 4 columns.
%      deltaq = is used when 'number_of_bins' or 'bin_width' is chosen. deltaq will represent
%               the number of bins when using 'number_of_bins' and will represent the bin 
%               width when using 'bin_width'.     
%  Function output = matrix of histograms M, built according to the sample and the chosen
%                    measuring operator: H_operator = 'center' or H_operator = 'integral'. 
%                    If H_operator = center, the output matrix M will be a three column 
%                    matrix (angle, measurement at the center of the bin, number of counts 
%                    in each bin). If H_operator = integral, the matrix M will be a four 
%                    column matrix (angle, left edge of the bin, right edge of the bin, 
%                    number of counts in each bin).

% Choosing default options. When only the two first input of M 
% are specified (numAngles and sample). In this case, we use the edges of the bins whose 
% width was chosen by Leonhardt's method and a measurement operator set to 'integral'.
if nargin == 2
  option = 'bin_leonhardt';
  H_operator = 'integral';
end

% The following line builds a vector of the same size as the output of the function 
% homodyne_samples. In https://arxiv.org/pdf/1805.07414.pdf this is 'N' (= 20,000).
num_measurements = size(samples, 1);

% In following line, we take the first column of samples to construct the structure 
% called QuadHist. All spaced angles are storaged in this structure.
angles = samples(1:numAngles);

% The quadrature measurement results, which can be the
% outputs of the function homodyne_samples, are stored in the matrix QuadHist(i).allQuads.
% [QuadHist(i).counts, QuadHist(i).edges] places these measurement results stored in 
% (QuadHist(i).allQuads) in bins. These bins are constructed using the histcounts function
% with the pre-determined bin width (in the four first cases of option) or the automatic
% binning algorithm (in the other cases of option).  
for i=1:numAngles
% Storing the measurement results in the matrix (QuadHist (i) .allQuads):
    QuadHist(i).angle = angles(i);
    QuadHist(i).allQuads = samples((i:numAngles:end),2);
% Building the histograms using a fixed number of bins, given by num_of_bins = deltaq:   
    if strcmp(option,'number_bins')
        num_of_bins = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, num_of_bins);
% Building the histograms using Scott's method to calculate the bin width:
    elseif strcmp(option,'scott_true')
        Bin_Width_Scott = 3.5*std(QuadHist(i).allQuads)*((num_measurements/numAngles)^(-1/3));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width_Scott);
% Building the histograms using a fixed bin width:
    elseif strcmp(option,'bin_width')
        Bin_Width = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
% Building the histograms using Leohnardt's equation to calculate the bin width:
    elseif strcmp(option,'bin_leonhardt')
        % n is an estimate of the mean photon number
        n = n_quadrature(samples,num_measurements);
        Bin_Width = pi/(2*sqrt(2*n+1));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
% Building the histograms using the automatic binning algorithm from the Matlab function 
% histcounts:
    elseif any(strcmp(option,{'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'}))
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinMethod', option);
     else 
         error('Error when the option is not (number_bins, bin_width,bin_leonhardt, scott_true, auto, scott, fd, integers, sturges or sqrt)')
    end
end

if strcmp(H_operator,'center')
    
    for i = 1:numAngles
        d = diff(QuadHist(i).edges)/2;
        QuadHist(i).centers = QuadHist(i).edges(1:end-1)+d;
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1), QuadHist(i).centers.',QuadHist(i).counts.'];
    end
   
elseif strcmp(H_operator,'integral')
    
    for i = 1:numAngles
        QuadHist(i).edgesArray = zeros(length(QuadHist(i).edges-1),2);
        QuadHist(i).edgesArray= [QuadHist(i).edges(1:end-1).',QuadHist(i).edges(2:end).'];
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1) , QuadHist(i).edgesArray ,QuadHist(i).counts.'];
    end
  
end

% The last step is to vertically concatenate the arrays in QuadHist.M:
M = vertcat(QuadHist.M);
M = M(M(:,end)> 0,:);

end