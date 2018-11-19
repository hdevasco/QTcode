function M = matrix_histogram(numAngles, samples, option, H_operator,deltaq)

% matrix_histogram discretiza os valores contínuos da saída da função homodyne_samples e retorna uma matriz com a discretização
%   Entradas da função:
%   numAngles = número de fases igualmente espaçadas de 0 a pi.
%   samples = matrix que contem os resultados das medições homódinas
%   option = indica o método para definir a largura da caixa do histograma
%           Se a opção for "number_of_bins", O código constrói o histograma
%            com o número de caixas indicadas por essa entrada. Se a opção for 
%            'bin_width', esse valor será exatamente a largura da caixa do histograma. 
%            Se a opção for 'bin_leonhardt', esse valor será exatamente a largura da caixa
%            da sugestão de Leonhardt. Caso contrário, as opções ('scott_true', 'auto', 'scott', 'fd', 'inteiros', 'sturges', 'sqrt')
%            indicarão a largura ideal da caixa do histograma de acordo com a distribuição de probabilidade na amostragem do estado quântico.

%   H_operator = indica a forma de como o operador de medição será escolhido.
%    Se "H_operator" for "center", o operador de medição a ser usado na tomografia representará a medição que ocorre no centro da caixa.
%    Se "H_operador" for "integral", o operador de medição a ser usado na tomografia representará a medida que ocorre ao longo do comprimento
%    de toda a caixa.

%   deltaq = indica a largura que você deseja escolher para a caixa do histograma.

% Saída da função = 
% Constrói a matriz de histogramas M de acordo com a amostra e o operador de medição que você deseja usar para a reconstrução do 
% estado: operador de medição representando a medição no centro da caixa (H_operator = center) ou medindo o operador ao longo do comprimento da 
% caixa (H_operator = integral). Se a opção escolhida para o uso do operador de medição for "center", M será uma matriz com 
% colunas (ângulo, medido no centro da caixa, número de contagens da caixa). Se a opção escolhida para o uso do operador de medição 
% for "integral", M será uma matriz com colunas (ângulo, borda esquerda da caixa, borda direita da caixa, número de contagens da caixa).

% A condição abaixo constrói a matriz do histograma usando as bordas das caixas cuja largura é escolhida pelo método Leonhardt usando o 
% operador de medição ao longo do comprimento da caixa, considerando que estes método para a escolha do operador nos dá uma melhor estimativa
% entre todos os casos analisados.

if nargin == 2
  option = 'bin_leonhardt';
  H_operator = 'integral';
end

% A linha abaixo constrói um vetor do tamanho da saída da função homodyne_samples
num_measurements = size(samples, 1);

% Utiliza a primeira linha da saída da função homodyne_samples para a construção da estrutura QuadHist.
angles = samples(1:numAngles);

% O próximo passo, constrói os histogramas a partir de uma matriz de estruturas usando as medidas homódinas da amostra(saída da função
% homodyne_samples) escolhendo o método para calcular a largura da caixa.

for i=1:numAngles
% constrói a matriz de estruturas a partir da amostra
    QuadHist(i).angle = angles(i);
    QuadHist(i).allQuads = samples((i:numAngles:end),2);
% a condição abaixo escolhe o número de caixas da discretização em cada fase  
    if strcmp(option,'number_bins')
        num_of_bins = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, num_of_bins);
% a condição abaixo utiliza o método de Scott para calcular a largura da caixa do histograma para a distribuição
% em cada fase
    elseif strcmp(option,'scott_true')
        Bin_Width_Scott = 3.5*std(QuadHist(i).allQuads)*((num_measurements/numAngles)^(-1/3));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width_Scott);
% a condição abaixo usa a largura desejada para calcular o comprimento da caixa do histograma
    elseif strcmp(option,'bin_width')
        Bin_Width = deltaq;
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
% a condição abaixo constrói a largura da caixa do histograma de acordo com a sugestão de Leohnardt
    elseif strcmp(option,'bin_leonhardt')
        n = n_quadrature(samples,num_measurements);
        Bin_Width = pi/(2*sqrt(2*n+1));
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinWidth', Bin_Width);
% a condição abaixo constrói a largura da caixa do histograma de acordo com a documentação histcounts do MATLAB
    elseif any(strcmp(option,{'auto', 'scott', 'fd', 'integers', 'sturges', 'sqrt'}))
        [QuadHist(i).counts, QuadHist(i).edges]=histcounts(QuadHist(i).allQuads, 'BinMethod', option);
     else 
         error('Error when the option is not (number_bins, bin_width,bin_leonhardt, scott_true, auto, scott, fd, integers, sturges or sqrt)')
    end
end

% Constrói M para que o operador de medição representa a medição que ocorre exatamente no centro da caixa do histograma.

if strcmp(H_operator,'center')
    
    for i = 1:numAngles
        d = diff(QuadHist(i).edges)/2;
        QuadHist(i).centers = QuadHist(i).edges(1:end-1)+d;
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1), QuadHist(i).centers.',QuadHist(i).counts.'];
    end
   
% Constrói M para que o operador de medição represente a medida que ocorre ao longo do comprimento da caixa do histograma.   

elseif strcmp(H_operator,'integral')
    
    for i = 1:numAngles
        QuadHist(i).edgesArray = zeros(length(QuadHist(i).edges-1),2);
        QuadHist(i).edgesArray= [QuadHist(i).edges(1:end-1).',QuadHist(i).edges(2:end).'];
        QuadHist(i).M = [repmat(angles(i),length(QuadHist(i).counts.'),1) , QuadHist(i).edgesArray ,QuadHist(i).counts.'];
    end
  
end
% O passo a seguir finaliza a execução do código, utiliza a matriz de estruturas para construir a matriz de histogramas M que possui a
% a discretização dos valores contínuos da amostra

M = vertcat(QuadHist.M);
M = M(M(:,end)> 0,:);

end
