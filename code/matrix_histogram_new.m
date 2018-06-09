function M = matrix_histogram_new(numAngles, samples, option, H_operator,deltaq)

method = {'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'};

num_measurements = size(samples, 1);

angles = samples(1:numAngles);

for i=1:numAngles
    QuadHist(i).angle = angles(i);
    QuadHist(i).allQuads = samples((i:numAngles:end),2);
    if option>8
        num_bins = option
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, num_bins);
    elseif option == 8
        Bin_Width_Scott = 3.5*std(QuadHist(i).allQuads)*((num_measurements/numAngles)^(-1/3));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width_Scott);
    elseif option == 7
        Bin_Width= deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
    elseif (option > 0) && (option < 7),
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinMethod', method{option});
    end
end

if strcmp(H_operator,'center'),
    
    for i = 1:numAngles;
        d = diff(QuadHist(i).edges)/2;
        QuadHist(i).centers = QuadHist(i).edges(1:end-1)+d;
        
     end
 end

if strcmp(H_operator,'integral')
    for i = 1:numAngles;
        QuadHist(i).edgesArray = zeros(length(QuadHist(i).edges-1),2);
        QuadHist(i).edgesArray= [QuadHist(i).edges(1:end-1).',QuadHist(i).edges(2:end).'];

    end
end

end
