%   Este código de tomografia é usado para simular medições
%   feito em um estado óptico de "cat-Schrodinger"para encontrar o máximo
%   Estado de verossimilhança para esse conjunto de medições.
%   Na reconstrução, usamos histogramas com medidas de quadratura das
%   respectivas fases de medição gerando uma matriz com linhas(ângulo,
%   centro do bin, número de contagens no bin) e a partir dessa matriz
%   fizemos um processo de optimização para calcular a fidelidade do estado
%   verdadeiro e o novo estado gerado pelo histograma.
%  Calculamos a fidelidade entre o estado estimado de máxima verossimilhança
%  com o estado puro (psi) assim como a fidelidade entre o estado estimado usando histograma  com o estado puro (psi).
%   O espaço de estado dimensional infinito para o oscilador harmônico será
%   Representado na base do número de fótons.
%   Vamos truncar o espaço de Hilbert em fótons maxPhotonNumber

% This tomography code is used to simulate measurements made in a "cat-Schrodinger" optical state
% to find the maximum likelihood state for this set of measurements.
% In the reconstruction, we use histograms with measures of quadrature of the respective measurement
% phases generating a matrix with lines (angle, center of the bin, number of counts in the bin) and from
% this matrix we made an optimization process to calculate the fidelity of the true state and The new state generated by the histogram.
% We calculated the fidelity between the estimated state of maximum likelihood with the pure state (psi) as well as the fidelity between
% the estimated state using histogram with the pure state (psi).
% The infinite dimensional state space for the harmonic oscillator will be represented
% on the basis of the number of photons.

tic;
clc;
clear;

%Maximum number of photons
Mph = [10];


% Values for the number of bins
b =[10, 50, 100, 200, 500, 800, 1000, 1250, 1500, 2000];

% Number of measurements
nM = [20000];

for i=1:length(Mph),
    for j=1:length(b),
        for k=1:length(nM),
            
            fileName = ['maxnumberMph',num2str(Mph(i)),'b',num2str(b(j)),'nM',num2str(nM(k)),'.mat'];
            
            if exist(fileName,'file') == 2,
                load(fileName);
            else
                maxPhotonNumber = Mph(i);
                numberbins = b(j);
                num_sim  = 100;
                F1 = zeros(num_sim,1);
                F2 = zeros(num_sim,1);
                F3 = zeros(num_sim,1);
                F4 = zeros(num_sim,1);
                W = zeros(num_sim,1);
                T1 = zeros(num_sim,1);
                T2 = zeros(num_sim,1);
                T3 = zeros(num_sim,1);
                T4 = zeros(num_sim,1);
                T5= zeros(num_sim,1);
               
                nMeasurements       = nM(k);
                etaDetector         = 0.9;
                maxIterations       = 2000;
                stoppingCriterion   = 0.01;
                alpha = 1;
                phase = 0;
                etaState = 0.8;
                t=1;
                
                save(fileName);
            end
            while t <= num_sim,
                
                fprintf('>> Generating angles... ');
                tic;
                
                m = 20;% Number of equally spaced angles
                angles = pi*[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,19]/m;
                angles = repmat(angles,1,ceil (nMeasurements/m));
                angles = angles(1:nMeasurements)';
                
                
                disp('>> init tables...');
                S = init_tables(maxPhotonNumber);
                
                disp('>> generate cat vector...');
                psi = generate_cat_vector(alpha, phase, S);
                
                disp('>> apply loss...');
                rho = apply_loss(psi,etaState,S);
                
                
                % Agora, fazemos as medições do estado rho. Observe que temos que especificar
                % Resultados de medição máximo e mínimo possíveis -7 e 7.
                
                % Now, we make the measurements of state rho.  Notice we have to specify
                % maximum and minimum possible measurement results -7 and 7.
                
                fprintf(['>> homodyne samples... iteraction ',num2str(t), '\n']);
                samples = homodyne_samples(-7,7,etaDetector,angles,rho,S);
                
                % Estrutura contendo o elemento POVM correspondente a cada medida
                % resultado. Observe que os POVMs nÃ£o sÃ£o projetores puros.
                % A eficiência do detector foi incluída no cálculo dos POVMs.
                
                % Structure containing the POVM element corresponding to each measurement
                % result.  Note that the POVMs are not pure projectors.  The homodyne
                % detector's efficiency has been included in the computation of the POVMs.
%                 tic;
                %     Povms = make_measurement_struct(samples,etaDetector,S);
                
                % Agora vamos usar o algoritmo R * rho * R até termos feito 2000 iteraÃ§Ãµes
                % Ou chegarmos a parar a diferença 0.01 (o que ocorrer primeiro).
                % A intercepção é um limite superior da diferença entre o verdadeiro
                % Log-verossimilhanÃ§a máxima e a log-verossimilhança do estado de iterações.
                
                % Now we will use the R*rho*R algorithm until we have done 2000 iterations
                % or we reach stoppingCriterion 0.01 (whichever happens first).
                % stoppingCriterion is an upper bound on the difference between the true
                % maximum log-likelihood and the log-likelihood of that iterations's state.
                
                
                
                
                %  [rhoML1, Diagnostics] = rrhor_optimization(samples, S, etaDetector, 0, maxIterations, stoppingCriterion, []);
                %  Saída é o estado de mÃ¡xima verossimilhança, rhoML e uma estrutura grande
                %   de diagnósticos contendo informações sobre o progresso de cada iteração.
                %
                %   Em vez de usar R * rho * R, podemos usar uma combinaÃ§Ã£o de R * rho * R seguido
                %   Por iterações de nosso novo algoritmo, a ascensão gradiente regularizada.
                %   Isso geralmente é mais rápido, especialmente se você quiser estar no
                %   Verdadeiro estado de máxima verossimilhança.
                
                % output is the maximum likelihood state, rhoML, and a big structure
                % Diagnostics containing information about each iterations's progress.
                
                % Rather than using R*rho*R, we can use a combination of R*rho*R followed
                % by iterations of our new algorithm, the regularized gradient ascent.
                % This is usually faster, especially if you want to be very close to the
                % true maximum likelihood state.
                fprintf(['>> combined optimization... iteraction ',num2str(t), '\n']);
                tic;
                [rhoML2, Diagnostics] = combined_optimization( samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                T1(t) = toc;
                % Fidelidade entre o estado verdadeiro e o estado estimado.
                
                % Fidelity of true state and estimate
%                 disp('>> fidelity (samples)...');
                tic;
                F1(t) = fidelity(rhoML2, rho);
                T2(t)= toc;
                
                
                %Constrói a matriz M que tem como linhas(Ângulo, centro do bin,nÚmero de contagens no bin)
                
                % It constructs the matrix M that has as lines (angle, center of the bin, number of counts in the bin)
                
                disp('>> matrix histogram ...');
                tic;
                M = matrix_histogram(samples,angles,nM,b(j),m);
                T3(t) = toc;
              
               
                %     Povmshistogram = make_measurement_struct(M,etaDetector,S);
                
                % [rhoML1, Diagnostics] = rrhor_optimization(M, S, etaDetector, 0, maxIterations, stoppingCriterion, []);
                
%                 fprintf(['>> combined optimization (M)... iteraction ',num2str(t), '\n']);
                tic;
                [rhohistogram, Diagnostics] = combined_optimization( M, S, etaDetector, 0, maxIterations, stoppingCriterion);
                T4(t) = toc;
                
                % Calcula a fidelidade entre o estado verdadeiro e o estado reconstruído
                % usando o histograma
                
                % Calculates the fidelity between the true state and the reconstructed state using the histogram
%                 disp('>> fidelity (rhohistogram x rho)...');
                 tic;

                F2(t) = fidelity(rhohistogram, rho);
                 T5(t)= toc;
                
                disp('>> fidelity (rhohistogram x psi)...');
                 
                  F3(t) = fidelity(rhohistogram, psi);
                
                
                disp('>> fidelity (rhoML2 x psi)...');
                
                F4(t) = fidelity(rhoML2,psi);
               
                
                W(t)= F1(t)-F2(t);
                w= mean(W);
                D = std(W);
                f1 = mean (F1); % Average fidelity (rhoML2,rho)
                f2 = mean (F2); % Average fidelity (rhohistogram, rho)
                f3 = mean (F3); % Average fidelity (rhohistogram, psi)
                f4 = mean (F4); % Average fidelity (rhoML2, psi)
                t1 = mean (T1); % Average time to build (rhoML2)
                t2 = mean (T2); % Average time to calculate  fidelity (rhoML2,rho)
                t3 = mean (T3); % Average time to build M
                t4 = mean (T4); %  Average time to build (rhohistogram)
                t5 = mean (T5); % Average time to calculate fidelity(rhohistogram, rho)
                d1 = std(T1); %  Standard deviation of the mean time of the rhoML2 construct
                d2 = std(T2); % Standard deviation of the fidelity estimation time (rhoML2, rho)
                d3 = std(T3); %  Standard deviation of the mean time of construction of M
                d4 = std(T4); % Standard deviation of the mean time of the rhohistogram construction
                d5 = std(T5); % Standard deviation of the fidelity estimation time (rhohistogram, rho)
                dF1= std(F1); % Standard deviation of the  fidelity of (rhoML2, rho)
                dF2= std(F2); % Standard deviation of the  fidelity of (rhohistogram, rho)
                home;
                fprintf('>> Progress: %.2f%%\n', t/num_sim*100);
                t=t+1;
                
                save(fileName);
            end
        end
    end
end
toc;
