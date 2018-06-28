function M = matrix_histogram_new(numAngles, samples, option, H_operator,deltaq)

num_measurements = size(samples, 1);

angles = samples(1:numAngles);

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
        Bin_Width= deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
    elseif any(strcmp(option,{'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'}))
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinMethod', option);
     else 
         error('Error when the option is not (number_bins, bin_width, scott_true, auto, scott, fd, integers, sturges or sqrt)')
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

M = vertcat(QuadHist.M);
M = M(M(:,end)> 0,:);

end
