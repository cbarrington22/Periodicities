% This USER-FRIENDLY SCRIPT contains code to test for periodicities in time series data using the method outlined in Section 2 of “Barrington et al. (2022) Wind speed as a dominant source of periodicities in reported emission rates of volcanic SO2”. 
 
% Created by C. BARRINGTON, October 2022.

% DESCRIPTION of script 
% Sections of the code which feature in "MasterScript.m" have been re-written to provide a user-friendly option for those wishing to apply similar analysis to their own time series data.
% Although it is not intended as a stand-alone script it includes lines of code for applying the (1) bootstrap approach, (2) Lomb-Scargle periodogram and for (3) plotting the PSD estimates between two variables. 
% This script differs slightly from the code used in “MasterScript.m” but the analysis procedure is the same.

%% HOW TO USE THIS SCRIPT 
% This script is divided into three sections which should be run in sequence: (1) Bootstrap resampling (2) Lomb-Scragle and (3) Plotting PSD estimate   
% Variables have been defined as default except for “timeseries” which must be defined before running any part of this code: 
% Define timeseries as: timeseries = [timestamp, parameter1, standard deviation1, parameter2, standard deviation2..]; 
timeseries="<define here>";
% If more than two parameters are included, “timeseries” may be expanded to include additional columns where the standard deviation is always included in next column over from the corresponding parameter.
% The timestamp must be continuous with constant sampling rate and included in the first column of “timeseries”. 
% The sampling rate “deltaT” is defined as the the number of datapoints the measuring instrument or protocol can provide per unit time. For the case of daily mean SO2 emission rate deltaT is 1 day. 
deltaT=1; 
% If the “parameters” are not uniformly sampled (or data is missing), the corresponding rows of "parameter" and "standard deviation" should contain “NaN”.
% Every “parameter” must have an accompanying “standard deviation” column. 
% Although variables may be re-defined, this code not been tested when the pre-defined variables differ from those used in Barrington et al. (2022). 
% All variables that should be tailored for the analysis are aligned to the left. Code which is indented is not intended to be changed. 


%% PART 1: BOOTSTRAP RESAMPLING
% We consider the daily mean SO2 emission rate which is reported, but to account for the large variation in SO2 emission rates recorded throughout the day we also incorporate the standard deviation 
% in daily SO2 emission rate in our approach We applied a bootstrapping approach to calculate the PSD from resampled time series of the daily SO2 emission rates 1000 times, each time selecting a 
% random value (mean +/- standard deviation). By doing so we assume the distribution of values is uniform for each day which is not always the case. 
% The following lines of code resample (by bootstrapping) a time series (timeseries), nRun number of times, between two values (l1) and (2) where l1 is the daily mean - standard deviation and l2 is the daily mean + standard deviation. 
 
% USER DEFINED VARIABLES 
nRuns=1000; % Set nRuns to the number of samples (typically at least 1000)
specCase=0; % Option for special case... see lines 51 to 60 
                               
 
                                % NOT INTENTED TO BE CHANGED
                                numPar=(width(timeseries)-1)/2; % Determines the number of parameters in the timeseries assuming the first column is the timestamp, and each parameter is accompanied by standard deviation in the subsequent column 
                                % Pre-allocates space to output variables 
                                bsMatrix=zeros(length(timeseries(:,1)),nRuns); % Creates a 2D output matrix equal to the length of the timeseries and number of samples to bootstrap 
                                bsMatrix(:,:,numPar)=zeros(length(timeseries(:,1)),nRuns); % Converts to 3D array, where the number of pages is equal to the number of parameters to bootstrap  
                                % Resamples timeseries data 
                                for par=1:numPar  
                                    ind=par*2; % Converts the parameter number to the column index  
                                    % Defines limits for resampling
                                    l1=timeseries(:,ind)-timeseries(:,ind+1); % Mean - standard deviation 
                                    l2=timeseries(:,ind)+timeseries(:,ind+1); % Mean + standard deviation 
                                    for bs = 1:nRuns % For all bootstrap samples 
                                        for dp=1:length(timeseries(:,1)) % For each datapoint in timeseries
                                            % Fill the datapoint (dp) value of the bootstrap sample (bs) for parameter (par) with random sample between M1 and M2
                                            bsMatrix(dp,bs,par)=l1(dp,:)+(l2(dp,:)-l1(dp,:)).*rand(1,1); 
                                            % If values are non-continuous, additional lines of code are needed here, for example in the case of plume direction 
                                            % ***SPECIAL CASE OF PLUME DIRECTION***
                                            % This parameter is not continuous (-180 to 180 deg.), i.e., plume direction of 176 deg. +/- 40 deg. IS NOT 216 to 136 but -144 to 136 
                                            % if par==specCase % If parameter is daily mean of used plume speed
                                                % if bsMatrix(dp,bs,par)<-180 % If the value is less than -180
                                                            % bsMatrix(dp,bs,par) = bsMatrix(dp,bs,par)+360; % Add 360 degrees to value 
                                                % elseif bsMatrix(dp,bs,par)>180 % If the value is greater than 180 
                                                             % bsMatrix(dp,bs,par) = bsMatrix(dp,bs,par)-360; % Subtract 360 degrees from value   
                                                % end 
                                            % end 
                                        end 
                                    end 
                                end 
 
                                
% Output is bsMatrix where its length is the length of the time series, the number of columns is now the number of samples where each page is one “parameter” 


%% PART 2: LOMB-SCARGLE  
% We use MATLAB’s ‘plomb’ function (MathWorks, 2021) to calculate the Lomb-Scargle PSD estimate of each bootstrapped time series, all of which have a sampling rate (deltaT) of one day. 
% We use an oversampling factor (ofac) of four, which resembles the zero-padding technique for FFT-based methods and is used to improve the resolution of the frequency grid over which the PSD is calculated. The maximum frequency is determined by 1/(2 x deltaT), which for uniformly spaced signals corresponds to the Nyquist frequency (Nyquist, 1928). We refer to this here 
% as the pseudo-Nyquist frequency, which represents the shortest period we are capable of detecting based on the deltaT (two-days). The minimum frequency is defined by 1/(ofac x N x deltaT), where 
% N is the number of days between the start and end dates of the time series and not necessarily the number of datapoints. We further limit the frequency when reporting the PSD estimates, 
% to consider periods between the pseudo-Nyquist frequency and N/4, corresponding to two-days and four-cycles per time series, e.g., 497-days at Mayon since the time series is 1991-days in length. 
% We therefore exclude periodicities with only two or three complete cycles. 
 
% USER DEFINED VARIABLES 
lim=4; % Defines low frequency limit – the minimum number of cycles in timeseries     
ofac=4; % Oversampling factor (4-8 may be reasonable)
% FALSE ALARM PROBABILITIES (FAPs)
Pfa=[50 10 1 0.01]/100; % False alarm probability levels of 50, 10, 1 and 0.01 (%)
% False alarm probability levels of 50, 10, 1 and 0.01 (%) correspond to significance levels of: 
% 50 % 
% 90 % 
% 99 %
% 99.99 % 
 
 
                                % NOT INTENTED TO BE CHANGED 
                                % Nyquist frequency (fn): fn = 1/2*deltaT 
                                Nf=1/(2*deltaT); % Nyquist frequency is 0.5 (i.e., one cycle per two days)
                                % 'plomb' function see: https://www.mathworks.com/help/signal/ref/plomb.html  
                                % Pxx is returned at round(fmax/fmin) frequency points where the minimum frequency is 1/(ofac x N x deltaT) - i.e., frequency grid is determined by the sampling interval and the length of the timeseries 
                                % N is the number of datapoints 
                                % fmax defaults to 1/(2*deltaT) which for uniformly spaced signals corresponds to Nf 
                                Pd=1-Pfa; % Returns the power-level threshold such that a peak with a value larger than output (pth) has a probability equal to Pd of being a true signal and not the result of a purely random signal
                                % Pre-allocates space to output variables
                                [~,ncol,npages]=size(bsMatrix);   % Determines number of columns and number of pages for 3D array (i.e. number of bootstrap runs and number of parameters)
                                % Runs Lomb-Scargle analysis once to obtain size of output parameters 
                                % For the first bootstrapped sample and the first parameter, calculates Lomb Scargle power spectral density and false alarm probability thresholds  
                                [pxx,f,pth]=plomb(bsMatrix(:,1,1),timeseries(:,1),Nf,ofac,'normalized','Pd',Pd); 
                                % Pre-allocates space for PSD output 
                                pxxMasterBs=zeros(length(pxx),ncol);pxxMasterBs(:,:,npages)=zeros(size(pxxMasterBs));
                                % Pre-allocates space for frequency grid output 
                                fMasterBs=zeros(length(f),ncol);fMasterBs(:,:,npages)=zeros(size(fMasterBs));
                                % Pre-allocates space for false alarm probability thresholds output 
                                pthMasterBs=zeros(length(pth),ncol);pthMasterBs(:,:,npages)=zeros(size(pthMasterBs));
                                % Pre-allocates space for calculated median and standard deviation 
                                msdMaster=zeros(length(pxx),2);msdMaster(:,:,npages)=zeros(size(msdMaster));
                                for par = 1:npages % For each parameter 
                                    for c = 1:ncol % And each bootstrapped sample
                                     % Calculates Lomb-Scargle periodogram and stores PSD, frequency grid and FAP threshold data
                                     % Set spectrum type to 'normalized' to get the standard Lomb-Scargle periodogram, which is scaled by two times the variance of x
                                        [pxxMasterBs(:,c,par),fMasterBs(:,c,par),pthMasterBs(:,c,par)]=plomb(bsMatrix(:,c,par),timeseries(:,1),Nf,ofac,'normalized','Pd',Pd); 
                                    end 
                                % For each parameter 
                                % Finds median and standard deviation for each frequency 
                                % Computes the median of each row  
                                msdMaster(:,1,par)=median(pxxMasterBs(:,:,par),2);
                                % Computes the standard deviation of each row 
                                msdMaster(:,2,par)=std(pxxMasterBs(:,:,par),0,2);
                                end 

                                % Checks returned frequency grid and FAP threshold are the same across parameters (sanity check) 
                                tf=isequal(fMasterBs(:,:,1),fMasterBs(:,:,2)); % Frequency grid 
                                tp=isequal(pthMasterBs(:,:,1),pthMasterBs(:,:,2)); % FAP thresholds
                                if tf==1&&tp==1 % If frequency grid and FAP thresholds are the same across parameters
                                % Saves variables as .mat file 
                                fprintf('Median and standard deviation of PSD calculated from bootstrapped parameters.\n'); 
                                else  
                                fprintf('Warning: Frequency grid and FAP thresholds for not equal for all parameters.');  
                                end
              
% Output is msdMaster where its length is the length of the time series, the two columns represent the median and standard deviation of the lomb-scargle applied to bootstrapped sampled, where each page is one “parameter” 


%% PART 3: PLOTTING PSD ESTIMATE
% We compare the PSD of SO2 emission rates and plume speeds for Mayon volcano between August 2011 and January 2017
% Refer to Fig. 6 and Fig. 7 in Barrington et al. (2022). 
                                      
% USER DEFINED VARIABLES 
par1=1; % Set “par1” to the y-axis parameter (Note: the first parameter in timeseries is parameter 1, the second will be parameter 2 etc.)
par2=2; % Set “par2” to the x-axis parameter
str1='<My first parameter>'; % y-axis parameter name 
str2='<My second parameter>'; % x-axis parameter name 
 
 
                                % NOT INTENTED TO BE CHANGED 
                                % Display settings 
                                dpCol=[0, 0.2, 0.9]; % Marker colour 
                                colSD=[0.5,0.5,0.5]; % Colour of error bars (standard deviation) 
                                maxP=length(timeseries)/lim;
                                medInd=1; % Column index for median
                                sdInd=2; % Column index for standard deviation
                                % Converts frequency data in fMasterBs to period 
                                t=ones(length(fMasterBs(:,1,1)),1); 
                                freq=fMasterBs(:,1,1); 
                                per=t./freq; % Converts to period
                                ofInt=find(per<=maxP); % Finds periods between the pseudo-Nyquist frequency and N/4
                                freq=per(ofInt,:); % Extracts periods between the pseudo-Nyquist frequency and N/4
                                % Figure  
                                figure;
                                e=errorbar(msdMaster(ofInt,medInd,par2),msdMaster(ofInt,medInd,par1),msdMaster(ofInt,sdInd,par1),msdMaster(ofInt,sdInd,par1),msdMaster(ofInt,sdInd,par2),msdMaster(ofInt,sdInd,par2),'.','CapSize',1);hold on; 
                                e.Color = colSD;
                                p = plot(msdMaster(ofInt,medInd,par2),msdMaster(ofInt,medInd,par1), '.');p.Color = dpCol;
                                % Adds FAP thresholds 
                                % Parameter 1
                                a=yline(pthMasterBs(1,1,par1),':');a.Color=[0.5 0 0.1];a.LineWidth=1;
                                b=yline(pthMasterBs(2,1,par1),'--');b.Color=[0.5 0 0.1];b.LineWidth=1;
                                c=yline(pthMasterBs(3,1,par1),'-');c.Color=[0.7 0 0];c.LineWidth=1;
                                d=yline(pthMasterBs(4,1,par1),'-');d.Color=[1 0 0];d.LineWidth=1;
                                % Parameter 2
                                a=xline(pthMasterBs(1,1,par2),':');a.Color=[0.5 0 0.1];a.LineWidth=1;
                                b=xline(pthMasterBs(2,1,par2),'--');b.Color=[0.5 0 0.1];b.LineWidth=1;
                                c=xline(pthMasterBs(3,1,par2),'-');c.Color=[0.7 0 0];c.LineWidth=1;
                                d=xline(pthMasterBs(4,1,par2),'-');d.Color=[1 0 0];d.LineWidth=1; 
                                % Set axis limits 
                                maxX=max(msdMaster(ofInt,medInd,par2));
                                maxY=max(msdMaster(ofInt,medInd,par1));
                                if maxX>pthMasterBs(1,1,1)||maxY>pthMasterBs(1,1,1) % If maximum of either axis is greater than the 50% threshold…
                                    if maxX>maxY % and x- axis is greater than y-axis.. 
                                        xlim([0 (maxX+max(msdMaster(ofInt,sdInd,par2)))+1]); % Set both axes to x-axis limit.. etc. 
                                        ylim([0 (maxX+max(msdMaster(ofInt,sdInd,par2)))+1]);  
                                    elseif maxY>maxX 
                                        xlim([0 (maxY+max(msdMaster(ofInt,sdInd,par1)))+1]); 
                                        ylim([0 (maxY+max(msdMaster(ofInt,sdInd,par1)))+1]); 
                                    end
                                    Fig = gca; cYLim = Fig.YLim; % Finds y-axis limit 
                                    for in = 1:length(pthMasterBs(:,1,1)) 
                                        if pthMasterBs(in,1,2)<max(cYLim)
                                           text(0.3*1,pthMasterBs(in,1,2)+.2,[num2str((((Pd(in))')*100)),repmat('%',[1,1])],'FontSize',11); % Adds FAP labels  
                                        end 
                                    end 
                                else % If max x- or y-axis isn't greater than the 50% threshold       
                                    xlim([0 pthMasterBs(1,1,1)+1]); % Set axis limit to show the 50% threshold
                                    ylim([0 pthMasterBs(1,1,1)+1]);  
                                    text(0.3* 1,pthMasterBs(1,1,2)+.2,[num2str((((Pd(1))')*100)),repmat('%',[1,1])],'FontSize',11);
                                    % end
                                end
                                xlabel(sprintf('Lomb-Scargle PSD estimate of %s',str1),'FontSize',11); % Axis labels according to user-defined string 
                                ylabel(sprintf('Lomb-Scargle PSD estimate of %s',str2),'FontSize',11);
                                Fig.FontSize=11;set(gcf,'color','w');

                                        
% Output is figure comparing the median PSD form two bootstrapped parameters.
 
% END OF SCRIPT    