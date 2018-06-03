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

% QuadHist = struct('numAngles', numAngles, ...
%     'QuadHist.angle', QuadHist.angle, ...
%     'QuadHist.allQuads', QuadHist.allQuads, ...
%     'QuadHist.counts', QuadHist.counts , ...
%     'QuadHist.edges',QuadHist.edges );

QuadHist = struct( numAngles, QuadHist.angle, QuadHist.allQuads, QuadHist.counts , QuadHist.edges);

if H_operator == 'center'
    for i = 1:numAngles
        M2 = zeros(1,3);
        d = diff(QuadHist(i).edges)/2;
        QuadHist(i).centers = QuadHist(i).edges(1:end-1)+d;
        
        if option > 8,
            
            C(:,i) = QuadHist(i).centers';
            A  = angles(i);
            N  = QuadHist(i).counts;
            A2 = A(:);
            C2 = C(:);
            N2 = N(:);
            M = [A2, C2, N2];
        else
            MA = [repmat(angles(i),length(QuadHist(i).counts),1), QuadHist(i).centers', QuadHist(i).counts'];
            M2 = [M2; MA];
        end
        M =M(2:end,:);
        ind = find(M(:,3) > 0);
        M = M(ind,:);
    end
end
if H_operator == 'integral'
    for i = 1:numAngles
        M2       = zeros(1,4);
        QuadHist(i).edgesArray = zeros(length(QuadHist(i).edges-1),2);
        QuadHist(i).edgesArray= [QuadHist(i).edges(1:end-1)',QuadHist(i).edges(2:end)'];
        MA = [repmat(QuadHist(i).angle ,length(QuadHist(i).counts(:,i)),1), QuadHist(i).edgesArray , QuadHist(i).counts(:,i)];
        M2 = [M2; MA];
    end
    M =M(2:end,:);
    ind = find(M(:,4) > 0);
    M = M(ind,:);
end

end
