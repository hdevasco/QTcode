function M = matrix_histogram(samples,angles,nM,b,m) 

% matrix_histogram constructs an array from samples whose columns are respectively
% . m equally spaced angles between 0 to pi with each angle being repeated a thousand times
% . Center of the boxes of each histogram formed by the quadrature measurements of samples of each angle separately, placed in "b" bins.
% . Number of quadrature measure counts in each bin in respective histograms

%Na matriz H, cada coluna é igual aos valores das medidas em quadratura da matriz de amostras para cada ângulo separadamente.
% In matrix H, each column equals the values of the quadrature measurements of the samples matrix for each angle separately.
H = zeros(nM/m,m);

% "A" será uma matriz coluna cujos ângulos repetidos equivalem as "m" fases igualmente espaçadas de 0 a pi.
% "A" will be a column matrix whose repeated angles are equal to "m" phases equally spaced from 0 to pi.

angles = pi*[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,19]/m;
angles = repmat(angles,b,(nM/m));
                angles = angles(1:b*m)';
            A = angles;
            
%Na matriz N, temos o número de contagem das caixas dos histogramas de cada ângulo igualmente espaçado.
% In matrix N, we have the count number of the bins of the histograms of each equally spaced angle.

N = zeros(b,m);

%Na matriz C, cada coluna é igual aos centros das caixas dos histogramas de cada ângulo separadamente, ambos colocados em 1000 caixas.
% In matrix C, each column equals the centers of the boxes of the histograms of each angle separately, both placed in 1000 bins.
C = zeros(b,m);

    for i=1:m,
         
        H(:,i) = samples((i:m:end),2); 
        
        [N(:,i),edges] = histcounts(H(:,i),b);
  
       d = diff(edges)/2;
       
  
       centers = edges(1:end-1)+d;
       
        C(:,i)= centers';
        
%         figure;
% 
%         title('h(i)');
%      
%         h(i) = histogram(H(:,i),b)
        
        
    end
    
    M = [A,C(:),N(:)];