tic;
clc;
clear;

maxPhotonNumber    = 10;
% maxPhotonNumber    = 20;

numMeasurements    = 20000;

option             = [10, 20, 40, 100];

std_fid_RhoML2          =         zeros(1,4);
mean_fid_RhoML2         =         zeros(1,4);
std_fid_RhoHistogram    =         zeros(1,4);
mean_fid_RhoHistogram   =         zeros(1,4);
std_fid_RhoScott        =         zeros(1,4);
mean_fid_RhoScott       =         zeros(1,4);

for c1=1:length(maxPhotonNumber),
    for c2=1:length(numMeasurements),
         for c3=1:length(option),
             
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(c1)),'op',num2str(option(c3)),'nM',num2str(numMeasurements(c2)),'.mat'];
  
            if exist(fileName1,'file') == 2,
                
                load(fileName1);  
            end 
             
            std_fid_RhoML2(c3)         = stdFML2;
            mean_fid_RhoML2(c3)        = meanFML2;
            std_fid_RhoHistogram(c3)   = stdFHistogram;
            mean_fid_RhoHistogram(c3)  = meanFHistogram;
            std_fid_RhoScott(c3)       = stdFScott;
            mean_fid_RhoScott(c3)      = meanFScott;

            
        end
    end
end
 
 hold on          
E1 = std_fid_RhoML2;
errorbar(option(1,4),mean_fid_RhoML2(1,4), E1(1,4),'k');

E2 = std_fid_RhoHistogram;
errorbar(option',mean_fid_RhoHistogram, E2,'r');

E3 = std_fid_RhoScott;
errorbar(option(1,4),mean_fid_RhoScott(1,4), E3(1,4),'b');

