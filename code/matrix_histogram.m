function M = matrix_histogram(numAngles, samples, option, H_operator,deltaq)

% matrix_histogram discretiza os valores contínuos da saída da função homodyne_samples e retorna uma matriz com a discretização
% Entradas da função:
% numAngles = número de fases igualmente espaçadas de 0 a pi.
% samples = matrix que contem os resultados das medições homódinas
% option = indica o método para definir a largura da caixa do histograma
%           Se a opção for "number_of_bins", O código constrói o histograma
%            com o número de caixas indicadas por essa entrada na função histcounts. Se a opção for 
%            'bin_width', esse valor será exatamente a largura da caixa do histograma a ser inserida em histcounts. 
%            Se a opção for 'bin_leonhardt', esse valor será exatamente a largura da caixa
%            da sugestão de Leonhardt, sendo essa implementada em histcounts. Se a opção for 'scott_true', a largura da caixa do 
%            histograma será escolhida usando a expressão analítica do método de Scott a ser implementada em histcounts. Caso contrário, as 
%            opções ('auto', 'scott', 'fd', 'inteiros', 'sturges', 'sqrt' indicarão a largura da caixa do histograma de acordo com 
%            a documentação histcounts do Matlab, construindo o histograma sem nenhuma largura previamente determinada. O tamanho do bin do 
%            histograma, nestes casos, buscam cobrir o intervalo de dados e revelar a forma da distribuição de probabilidade na amostragem 
%            do estado quântico.

%   H_operator = indica a forma de como o operador de medição será escolhido.
%    Se "H_operator" for "center", o operador de medição a ser usado na tomografia representará a medição que ocorre no centro da caixa.
%    Se "H_operador" for "integral", o operador de medição a ser usado na tomografia representará a medida que ocorre ao longo do comprimento
%    de toda a caixa.

%   deltaq = poderá indicar o número de caixas ou a largura de entrada que se deseja construir a caixa do histograma. Se 'option' for 
%   "number of bins", deltaq representará o número de caixas do histograma. Se 'option' for "Bin_Width", deltaq representará a largura que se 
%    deseja ter para a caixa do histograma.

% Saída da função = 
% Constrói a matriz de histogramas M de acordo com a amostra e o operador de medição que você deseja usar para a reconstrução do 
% estado: operador de medição representando a medição no centro da caixa (H_operator = center) ou medindo o operador ao longo do comprimento da 
% caixa (H_operator = integral). Se a opção escolhida para o uso do operador de medição for "center", M será uma matriz com 
% colunas (ângulo, medida no centro da caixa, número de valores de medição presentes em cada caixa). Se a opção escolhida para o uso do 
% operador de medição for "integral", M será uma matriz com colunas (ângulo, borda esquerda da caixa, borda direita da caixa, 
% número de valores de medição presentes em cada caixa).

% A condição abaixo constrói a matriz do histograma usando as bordas das caixas cuja largura é escolhida pelo método Leonhardt usando o 
% operador de medição ao longo do comprimento da caixa, considerando que estes método para a escolha do operador nos dá uma melhor estimativa
% entre todos os casos analisados.

if nargin == 2
  option = 'bin_leonhardt';
  H_operator = 'integral';
end

% A linha abaixo constrói um vetor do tamanho da saída da função homodyne_samples
num_measurements = size(samples, 1);

% A linha 50 utiliza a primeira coluna de samples para a construção da estrutura QuadHist, todos os ângulos igualmente espaçados são 
% armazenados nessa estrutura.
angles = samples(1:numAngles);

% Da linha 62 a linha 79, as linhas que possuem ([QuadHist(i).counts, QuadHist(i).edges]) constroem os histogramas usando histcounts a partir
% da matriz de estruturas (QuadHist(i).allQuads), essa estrutura guarda todos os valores de medida de quadratura de samples(saída da função
% homodyne_samples).

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
