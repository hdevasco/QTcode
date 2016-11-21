function h = homodyne_hist(samples,bidWidth)
X = randn(samples,1);% Distribute the samples in the bins
% edges = [lim inf x x x x ...x lim sup]; The first element is the left edge of the first share (bin 1) and the last element is the right edge of the last share (bin "n") of the border vector.
N = histcounts(X,edges)
end