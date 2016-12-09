function M = matrix_histogram(samples) 
% matrix_histogram constructs an array from samples whose columns are respectively
% . m equally spaced angles between 0 to pi with each angle being repeated a thousand times
% . Center of the boxes of each histogram formed by the quadrature measurements of samples of each angle separately, placed in "b" bins.
% . Number of quadrature measure counts in each bin in respective histograms

    for (i=1:m)
        
        A(:,i) = samples((i:m:end),1);
        
        H(:,i) = samples((i:m:end),2); 
        
         [N(:,i),edges] = histcounts(H(:,i),b);
  
       d = diff(edges)/2;
       
  
       centers = edges(1:end-1)+d;
       
        C(:,i)= centers';
        
        figure;

        title('h(i)');
     
        h(i) = histogram(H(:,i),b)
        
    end