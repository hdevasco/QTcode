function M = matrix_histogram(samples,b,m) 
% matrix_histogram constructs an array from samples whose columns are respectively
% . m equally spaced angles between 0 to pi with each angle being repeated a thousand times
% . Center of the boxes of each histogram formed by the quadrature measurements of samples of each angle separately, placed in "b" bins.
% . Number of quadrature measure counts in each bin in respective histograms

%Na matriz H, cada coluna é igual aos valores das medidas em quadratura da matriz de amostras para cada ângulo separadamente.
% In matrix H, each column equals the values of the quadrature measurements of the samples matrix for each angle separately.
H = zeros(b,m);

%Na matriz A, cada ângulo é repetido "b" vezes em cada coluna para transformá-la em um vetorcoluna para a construção de M.
% In matrix A, each angle is repeated a thousand times in each column in order to transform it into a column vector for the construction of M.
A = zeros(b,m);

%Na matriz N, temos o número de contagem das caixas dos histogramas de cada ângulo igualmente espaçado.
% In matrix N, we have the count number of the bins of the histograms of each equally spaced angle.

N = zeros(b,m);

%Na matriz C, cada coluna é igual aos centros das caixas dos histogramas de cada ângulo separadamente, ambos colocados em "b" caixas.
% In matrix C, each column equals the centers of the boxes of the histograms of each angle separately, both placed in "b" bins.
C = zeros(b,m);

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
  
  M=[A(:),C(:),N(:)];