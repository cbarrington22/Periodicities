% This MASTER SCRIPT has been used to produce “Barrington et al. (2022) Wind speed as a dominant source of periodicities in reported emission rates of volcanic SO2” and calls on the functions contained in “MasterScript_functions.zip” to analyse time series data from the Network of Observation of Volcanic and Atmospheric Change (NOVAC) database for the purpose of identifying periodic trends in daily observed SO2 emission rates.

% Created by C. BARRINGTON, March 2022.
% Modifications have not been logged. 

% DISCLAIMER: This MASTER SCRIPT has been provided to in the interest of research integrity. 
% This version of the code is not thoroughly commented nor is it free of redundancies and figures produced using this code have been modified and formatted independently although the data remains unchanged. 
% Some lines of code (including those within functions) are required to be commented/uncommented depending on whether the desired figure contains results from a single volcano or multiple subplots. It may therefore be necessary to run each section of the code separately. 
% Not all hardcoding has been removed. 
% Users wishing to apply similar analysis to their own datasets, are instead directed to the USER-FRIENDLY SCRIPT ("UserFriendly.m") which provides the necessary code to apply the method to other time series data. 

% DESCRIPTION of script 
% This script uses the following DATA:
 
% Data files from NOVAC database (https://novac.chalmers.se/):
% Column 1 - Date_UT_[YYMMDD]:                              Date of measurements
% Column 2 - SO2_flux_(daily_mean)_[kg s-1]:                Daily mean of observed SO2 flux
% Column 3 - SO2_flux_(daily_s.d.)_[kg s-1]:                Daily s.d of observed SO2 flux
% Column 4 - SO2_flux_(daily_25%)_[kg s-1]:                 Daily first quartile of observed SO2 flux
% Column 5 - SO2_flux_(daily_50%)_[kg s-1]:                 Daily second quartile of observed SO2 flux
% Column 6 - SO2_flux_(daily_75%)_[kg s-1]:                 Daily third quartile of observed SO2 flux
% Column 7 - SO2_plume_speed_(daily_mean)_[m s-1]:          Daily mean of used plume speed
% Column 8 - SO2_plume_speed_(daily_s.d.)_[m s-1]:          Daily s.d. of used plume speed
% Column 9 - SO2_plume_direction_(daily_mean)_[deg]:        Daily mean of observed plume direction
% Column 10 - SO2_plume_direction_(daily_s.d.)_[deg]:       Daily s.d. of observed plume direction
% Column 11 - SO2_plume_altitude_(daily_mean)_[m asl]:      Daily mean of observed plume altitude
% Column 12 - SO2_plume_altitude_(daily_s.d.)_[m asl]:      Daily s.d. of observed plume altitude
% Column 13 - SO2_plume_distance_(daily_mean)_[m]:          Daily mean of observed distance to plume
% Column 14 - SO2_plume_distance_(daily_s.d.)_[m]:          Daily s.d. of observed distance to plume
% Column 15 - SO2_plume_width_(daily_mean)_[m]:             Daily mean of plume width scanned by the instrument(s)
% Column 16 - SO2_plume_width_(daily_s.d.) [m]:             Daily s.d. of plume width scanned by the instrument(s)
% Column 17 - Cloud_cover_(daily_mean)_[%]:                 Daily mean of modeled cloud cover at summit altitude
% Column 18 - Cloud_cover_(daily_s.d.)_[%]:                 Daily s.d. of modeled cloud cover at summit altitude
% Column 19 - Total_number_of_measurements:                 Number of total scan measurements on each day
% Column 20 - Valid_number_of_measurements:                 Number of valid SO2 flux measurements on each day    
 
% Eruption history from GVP (https://volcano.si.edu/):
% Column 1 - Start Date
% Column 2 - Stop Date
% Column 3 - Eruption Certainty
% Column 4 - VEI
% Column 5 - Evidence
% Column 6 - Activity Area or Unit
 
% Script BREAKDOWN: 
% Runs 'plotTimeseries' function to load data files and if more than 8-days of data exist (based on minimum period 2-day and minimum number of cycles 4), PLOTS TIME SERIES OF SO2 EMISSIONS (daily mean +/- standard deviation) with volcanic activity 
        % Saves the following information <volcano>.mat in the directory 'outDirMat': 
            % Start date, end date and number of datapoints included in NOVAC database 
            % Dataset including: 
                % 1. Date 
                % 2. Daily mean of observed SO2 flux
                % 3. Daily s.d of observed SO2 flux
                % 4. Daily first quartile of observed SO2 flux
                % 5. Daily second quartile of observed SO2 flux
                % 6. Daily second quartile of observed SO2 flux
                % 7. Daily mean of used plume speed
                % 8. Daily s.d. of used plume speed
                % 9. Daily mean of observed plume direction
                % 10. Daily s.d. of observed plume direction
                % 11. Daily mean of observed plume altitude
                % 12. Daily s.d. of observed plume altitude
                % 13. Daily mean of observed distance to plume
                % 14. Daily s.d. of observed distance to plume
                % 15. Daily mean of plume width scanned by the instrument(s)
                % 16. Daily s.d. of plume width scanned by the instrument(s)
                % 17. Daily mean of modeled cloud cover at summit altitude
                % 18. Daily s.d. of modeled cloud cover at summit altitude
                % 19. Number of total scan measurements on each day
                % 20. Number of valid SO2 flux measurements on each day
         % Function saves figure in the directory 'outDirFig', <volcano_timeseriesSO2.fig>
    
% Runs 'plotSynthetic' function to CREATES SERIES OF SYNTHETIC SIGNALS with known periodicities which are resampled at the same time-points as the SO2 emission data at the selected volcano(es)
        % Four synthetic signals are created each including noise (using snr) and resampled to the range of SO2 emissions observed at the selected volcano(es):
        % Syn 1: with known periodicity equal to input value p1
        % Syn 2: with known periodicities equal to input values p2 and p3
        % Syn 3: with known periodicity equal to input values p1, p2 and p3
        % Random: no periodicity 
        % Function saves figure of each synthetic signal in the directory 'outDirSyn', <volcano_synthetic<number>.fig>
        
% Runs 'lombScargleSyn' function to CALCULATE LOMB-SCARGLE POWER SPECTRAL DENSITY (PSD) for each synthetic signal and false alarm probabilities defined by Pd
        % This function uses the 'plomb' function with the maximum frequency (Nf) and oversampling factor (ofac) defined 
        % This function saves the Lomb-Scargle periodogram for each synthetic signal in the directory 'outDirFig', <volcano_synthetic<number>_LS.fig>
 
% Runs 'plotLSSyn' function to PLOT PSD FOR SYNTHETIC SIGNALS  
        % Converts frequency data to period (days)
        % Displays periodicities within Nf and lim (for example, 2 days to number of days equal to 4 cycles per time series) 
        % Function saves figure in the directory 'outDirSyn', <volcano_SynPsdLS.fig>  
 
% Runs 'bootStrap' function to loop through all parameters in the dataset which are reported as a daily mean with standard deviation and RESAMPLES them using random numbers between the minimum (daily mean - standard deviation) and maximum (daily mean + standard deviation) values  
        % Function saves resampled data as 3D matrix, <volcano_bsMatrix.mat>:
            % Rows: each day (length of time series)
            % Column: each run (defined by nRuns)
            % Page: each parameter where:
                % 1: Daily mean of observed SO2 flux
                % 2: Daily mean of used plume speed
                % 3: Daily mean of observed plume direction
                % 4: Daily mean of observed plume altitude
                % 5: Daily mean of observed distance to plume
                % 6: Daily mean of plume width scanned by the instrument
                % 7: Daily mean of modeled cloud cover at summit altitude
            % This function saves non-resampled data as a 2D matrix, <volcano_nbsMatrix.mat>:
                % Number of total scan measurements on each day
                % Number of valid SO2 flux measurements on each day
                
% Runs 'lombScargleBs' function to CALCULATES LOMB-SCARGLE POWER SPECTRAL DENSITY OF BOOTSTRAPPED DATA (PSD) and FALSE ALARM PROBABILITIES (FAP) defined by Pd for all parameters in <volcano_bsMatrix.mat> 
       % This function uses the 'plomb' function with the maximum frequency (Nf) and oversampling factor (ofac) defined 
            % It checks the returned frequency grid and false alarm probability thresholds are equal across parameters 
            % This function calculates the PSD for all parameters contained in <volcano_bsMatrix.mat> and saves it in <volcano_psdLS.mat> 
            % This function calculates the median and standard deviation at each frequency for parameters contained in <volcano_bsMatrix.mat> and saves the output in <volcano_msdLS.mat> 
        % Note: A warning has been suppressed which indicates a parameter has no variance. This is likely due to a constant value for the plume altitude being used but may not be the case for all volcanoes
 
% Runs 'lombScargleNbs'function to CALCULATES LOMB-SCARGLE POWER SPECTRAL DENSITY OF NON-BOOTSTRAPPED DATA (PSD) and FALSE ALARM PROBABILITIES (FAP) defined by Pd for all parameters in <volcano_nbsMatrix.mat>
       % This function uses the 'plomb' function with the maximum frequency (Nf) and oversampling factor (ofac) defined 
            % It checks the returned frequency grid and false alarm probability thresholds are equal across parameters 
            % This function calculates the PSD for all parameters contained in <volcano_nbsMatrix.mat> and saves it in <volcano_psdLS.mat>         
        
% Runs 'plotLsSO2' function PLOTS MEDIAN PSD for SO2 EMISSION WITH STANDARD DEVIATION 
        % This function saves a variable containing the frequencies at which the PSD of SO2 emission crosses the corresponding Pd value (sigV)
        % Converts frequency data to period (days)
        % Displays periodicities within Nf and lim (for example, 2 days to number of days equal to 4 cycles per time series) 
        % Function saves figure in the directory 'outDirFig', <volcano_LsSO2.fig>
        
% Runs 'plotLsAllP' function PLOTS MEDIAN PSD WITH STANDARD DEVIATION FOR ALL PARAMETERS included in <volcano_msdLS.mat> and calculated PSD for parameters included in <volcano_psdLS.mat>
        % This function plots an x-line at the frequencies contained in sigV, across all plots 
        % Function saves figure in the directory 'outDirFig', <volcano_LsAll.fig>
        
% Runs 'plotPSDplumespeed' function to PLOTS THE PSD ESTIMATES OF SO2 EMISSION RATE AGAINST THE PSD ESTIMATES OF PLUME SPEED
        % Function saves figure in the directory 'outDirFig', <volcano_PSDxy.fig>
    
        
             % HARDCODING 
                    % WORKSPACES, DIRECTORIES AND FILENAMES 
                    % Adds path to project WORKSPACE 
                    addpath <insert path to project workspace> % Project workspace
                    addpath <insert path to functions> % Access to functions 
                    addpath <insert path to data files>% Path to datafiles
                    addpath <insert path to created .mat files> % Access to created .mat files 
 
                    % Turns off WARNINGS 
                    warning('off', 'MATLAB: datenum:EmptyDate');
                    warning('off','signal:plomb:SignalZeroVariance'); % Warning indicating 'zero variance' in time series (occurs when plume altitude is unknown and assumed to be summit altitude) 
 
                    % Sets DIRECTORIES 
                    dirNovac = '<insert path to local NOVAC data>'; % Directory to .csv files from NOVAC database
                    dirGvp = '<insert path to local GVP data>'; % Directory to eruptive history modified from GVP
                    % GVP data: 
                        % Eruptive history was downloaded from GVP for each volcano 
                        % Data before Jan 1st 2000 has been excluded 
                        % Unconfirmed eruptions have been removed 
                        % Where no eruption end date is listed, the eruption start date has been used, unless ongoing (i.e., latest entry in database) 
                    outDirMat = '<insert path to created .mat files >'; % Directory to output folder for saving created data files
                    outDir = '<insert path to output directory>‚Äô; % Directory to output folder for saving results file';
                    outDirSyn = '<insert path to output directory for synthetic data>'; % Directory to output folder for saving datafiles and results from synthetic data 
                    outDirFig = '<insert path to main figure directory>'; % Directory to output folder for saving figures
                    outDirFigBS = '< insert path to un-used figure directory>'; % Directory to output folder for saving resampling figures (sanity check only)
 
                    % FILENAMES 
                    inFile = '<insert name of .mat file containing structure with volcano names >'; % Directory to load structure which lists the name of the volcanoes to analyse
                    ext1bs = '_bsMatrix'; % Filename extension for bootstrapped data 
                    ext2bs = '_msdLS'; % Filename extension for Lomb-Scargle results of bootstrapped data 
                    ext1nbs = '_nbsMatrix'; % Filename extension for non-bootstrapped data 
                    ext2nbs = '_psdLS'; % Filename extension for Lomb-Scargle results of non-bootstrapped data 
                    extSyn = '_SynPsdLS'; % Filename extension for Lomb-Scargle results of synthetic data 
 
                    % FILE TYPES
                    fileTypeNovac = '.txt'; % File type for local NOVAC data file 
                    fileTypeGvp = '.csv'; % File type for local GVP data file 
                    fileTypeMat = '.mat'; % File type for data files created in this script 
                    % DATE FORMATS
                    formatInNovac = 'yyyymmdd'; % Defines date format used in NOVAC 
                    formatInGVP = 'yyyymmmdd'; % Defines date format used in GVP data 
 
            % PRE-DEFINED VARIABLES 
                    % ERUPTIVE HISTORY 
                    minWidth = 10; % Minimum number of days used to indicate eruption if single event
 
                    % SYNTEHTIC signals 
                    % Periodicities (days)
                    p1 = 120; % ~4-months periodicity
                    p2 = 50; % 50-day periodicity
                    p3 = 14; % 14-days periodicity
                    % 'awgn' function see: https://www.mathworks.com/help/comm/ref/awgn.html
                    snr = 0.01; % Signal to noise ratio 10% of signal 
 
                    % BOOTSTRAPPING 
                    nRuns = 1000; % Number of times time series is resampled
 
                    % LOMB-SCARGLE  
                    deltaT = 1; % Sampling rate is one measurement per day 
                    % Nyquist frequency (fn): fn = 1/2*deltaT 
                    Nf = 1/(2*deltaT); % Nyquist frequency is 0.5 (i.e., one cycle per two days)
                    lim = 4; % Defines low frequency limit, minimum number of cycles in timeseries     
                    % 'plomb' function see: https://www.mathworks.com/help/signal/ref/plomb.html 
                    % Pxx is returned at round(fmax/fmin) frequency points where the minimum frequency is 1/(ofac x N x deltaT) - i.e., frequency grid is determined by the sampling interval and the length of the timeseries 
                    % N is the number of datapoints 
                    % fmax defaults to 1/(2*deltaT) which for uniformly spaced signals corresponds to Nf 
                    % To improve the resolution, more frequencies can be tested by incorporating an oversampling factor, ofac (4-8 is reasonable), default is 4
                    % The use of ofac to interpolate or smooth a spectrum resembles the zero-padding technique for FFT-based methods
                    ofac = 4; % Oversampling factor
 
                    % FALSE ALARM PROBABILITIES (FAPs)
                    Pfa = [50 10 1 0.01]/100; % False alarm probability levels of 50, 10, 1 and 0.01 (%)
                    Pd = 1 - Pfa; % Returns the power-level threshold such that a peak with a value larger than output (pth) has a probability equal to Pd of being a true signal and not the result of a purely random signal
                    sigV = 3; % Determines the significance level in Pd to determine the plotting of x-line. For example, sigV = 3, will plot an x-line at periodicities in SO2 emission where displayed PSD crosses pd(:,3)
                    % False alarm probability levels of 50, 10, 1 and 0.01 (%) correspond to significance levels of: 
                    % 50 % 
                    % 90 % 
                    % 99 %
                    % 99.99 % 
 
                    % DISPLAY
                    % COLOURS for plotting
                    colSO2 = [0, 0.2, 0.9]; % SO2
                    colplumespeed = [0, 0.7, 0]; % Plume speed
                    colplumedir = [0, 0.5, 0]; % Plume direction 
                    colalt = [0, 0.1, 0.3]; % Plume altitude 
                    coldist = [0.6, 0.2, 0.5]; % Distance to plume
                    colwidth = [0.5,0.4,1]; % Plume width 
                    colcloud = [0.4, 0.4, 0.5]; % Cloud cover
                    coltotalnum = [1, 0.8, 0]; % Total measurements 
                    colvalmeas = [0.9, 0.5, 0.1]; % Valid measurements 
                    colVEI = [1, 0, 0, 0.5]; % VEI 
                    colSyn = [0.7, 0, 0.5]; % Synthetic signal 
                    colSD = [0.5,0.5,0.5]; % Standard deviation (PSD) 
                    colSel = [colSO2; colplumespeed; colplumedir; colalt; coldist; colwidth; colcloud; coltotalnum; colvalmeas]; % Compiled for use in function
 
                    % Figure SIZE 
                    x0 = 10;
                    y0 = 10;
                    wX = 1000;
                    hY = 500;
     
 
% EVALUATION ROUTINE 
% Loads list of volcanoes to analyse 
novacDatabase = ([dirNovac, inFile]); novacDatabase = load(novacDatabase); novacDatabase = novacDatabase.novacDatabase; % Loads structure listing names of volcanoes to analyse 
volcano = fieldnames(novacDatabase); % Creates variable with volcano names 
 
% Volcano index: 
% 1 Arenal
% 2 Concepcion
% 3 Copahue
% 4 Cotopaxi
% 5 Etna
% 6 FuegoDeColima
% 7 FuegoGuatemala
% 8 Galeras
% 9 Isluga
% 10 Lascar
% 11 Llaima
% 12 Masaya
% 13 Mayon
% 14 Momotombo
% 15 NevadoDelRuiz
% 16 Nyiragongo
% 17 PitonDeLaFournaise
% 18 PlanchonPeteroa
% 19 Popocatepetl
% 20 Sabancaya
% 21 SanCristobal
% 22 Sangay
% 23 SanMiguel
% 24 SantaAna
% 25 Santiaguito
% 26 Sinabung
% 27 Telica
% 28 Tungurahua
% 29 Turrialba
% 30 Ubinas
% 31 Villarrica
% 32 Vulcano
 
% TIME SERIES OF SO2 EMISSION RATE 
% plotTimeseries
% Loads NOVAC datafile and creates time series with all days between start and end date
% Plots daily SO2 emission rate with volcanic activity 
asp = 1; figure % Comment or not depending on whether plotting figure for one volcano or several (function might need amending)
for volIndx = 1:length(volcano) % For all volcanoes in list 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrNovac =[vName fileTypeNovac]; inDirNovac = fullfile(dirNovac,fileNameStrNovac); % Determines full directory to access data 
    [asp] = plotTimeseries(inDirNovac, vName, formatInNovac, fileTypeGvp, dirGvp, formatInGVP, colSO2, minWidth, colVEI, x0, y0, wX, hY, outDirFig, outDirMat, asp);
end                         
fprintf('Datafiles created!\n'); 
fprintf('Working on synthetic..\n'); 
 
% SYNTHETIC
for volIndx = [2, 13] % For Concepcion and Mayon  
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat =[vName fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2
        load(fileNameStrMat); 
        [sOut] = plotSynthetic(vName, dataSet, p1, p2, p3, snr, colSO2, colSyn, outDirSyn, outDirMat);
        % Note: check returned frequency grid and FAP thresholds
        [pxxMaster, fMaster, pthMaster] = lombScargleSyn(sOut,dataSet, Nf, ofac, Pd, vn, outDirSyn, extSyn);
        plotLSSyn(pxxMaster, fMaster, pthMaster, dataSet, deltaT, lim, p1, p2, p3, snr, colSyn, Nf, Pd, vn, outDirSyn, extSyn,  x0, y0, hY, wX);
    else
        fprintf('There is no corresponding .mat file for %s\n', vName); 
    end 
    fprintf('Synthetic complete!\n'); 
end
 
% BOOTSTRAPPING
for volIndx = 1:length(volcano) % For all volcanoes in list 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat = [vName fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2
        load(fileNameStrMat); 
        bootStrap(vn, dataSet, nRuns, outDirFigBS, outDirMat, ext1bs, ext1nbs);
    else
        fprintf('There is no corresponding .mat file for %s\n', vName); 
    end 
end 
        
% LOMB-SCARGLE 
for volIndx = 1:length(volcano) % For all volcanoes 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    % Bootstrapped parameters
    fileNameStrMat = [vName ext1bs fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2
        load(fileNameStrMat); % Using bsMatrix
        % Note: check on returned frequency grid and FAP thresholds
        lombScargleBs(bsMatrix, dataSet, Nf, ofac, Pd, vn, outDirMat, ext2bs)
    else
        fprintf('Skipping %s\n', vName); 
    end 
    fprintf('Lomb-Scargle complete for bootstrapped parameters! Working on non-bootstrapped parameters..\n') 
    % Non-bootstrapped samples 
    fileNameStrMat = [vName ext1nbs fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2
        load(fileNameStrMat); % Using nbsMatrix
        % Note: check on returned frequency grid and FAP thresholds
        lombScargleNbs(nbsMatrix, dataSet, Nf, ofac, Pd, vn, outDirMat, ext2nbs)
    else
        fprintf('Skipping %s\n', vName); 
    end 
end
fprintf('Lomb-Scargle complete for all volcanoes!\n') 
 
% PLOTS PSD (SO2)
fap = [];
asp = 1; figure % Comment or not depending on whether plotting figure for one volcano or several (function might need amending)
for volIndx = 13 % 1:length(volcano) % For Mayon only OR; % For all volcanoes 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat = [vName ext2bs fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2
        load(fileNameStrMat); 
        [asp, fap] = plotLsSO2(vn, fMasterBs, dataSet, msdMaster, colSD, colSO2, pthMasterBs, lim, deltaT, Nf, Pd, outDirFig, x0, y0, 500, hY, asp, fap); % SO2
    else
        fprintf('LS Periodogram not plotted for %s SO_2 emission rate\n', vName); 
    end 
end
 
% PLOTS PSD (ALL PARAMETERS)
for volIndx = 13 % 1:length(volcano) % For Mayon only OR; % For all volcanoes 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat1 = [vName ext2bs fileTypeMat]; exOut1 = exist(fileNameStrMat1, 'file'); 
    fileNameStrMat2 = [vName ext2nbs fileTypeMat]; exOut2 = exist(fileNameStrMat2, 'file');
    if exOut1 == 2 && exOut2 == 2
        load(fileNameStrMat1, 'dataSet', 'fMasterBs', 'pthMasterBs', 'msdMaster'); % Bootstrapped parameters
        load(fileNameStrMat2, 'pxxMasterNbs', 'fMasterNbs', 'pthMasterNbs'); % Non-bootstrapped parameters
        plotLsAllP(vn, fMasterBs, dataSet, fMasterNbs, pthMasterNbs, pthMasterBs, msdMaster, pxxMasterNbs, colSel, deltaT, Nf, colSD, lim, x0, y0, wX, hY, outDirFig);
    else
        fprintf('LS Periodogram not ploted for %s all parameters\n', vName); 
    end 
end
 
% PLOTS PSD OF SO2 EMISSION RATE AND PLUME SPEED
asp = 1; figure
for volIndx = 18:length(volcano) % [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17]  %1:length(volcano) % For all volcanoes 
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat = [vName ext2bs fileTypeMat]; exOut = exist(fileNameStrMat, 'file'); 
    if exOut == 2 
        load(fileNameStrMat, 'dataSet', 'fMasterBs', 'pthMasterBs', 'msdMaster'); % Bootstrapped parameters
        [asp] = plotPSDplumespeed(msdMaster, pthMasterBs, colSD, colSO2, Pd, vn, outDirFig, dataSet, deltaT, Nf, lim, fMasterBs, asp);
    else
        fprintf('PSD plume speed against PSD SO_2 emission rate not ploted for %s all parameters\n', vName); 
    end 
end
fprintf('Analysis complete!\n'); 
 
% APPENDIX D: PLOTS PDS FOR ALL PARAMETERS
for volIndx = [5, 6, 8, 12, 15, 16, 21, 28, 29, 31] 
    figure 
    subplot(2,1,1)
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrNovac =[vName fileTypeNovac]; inDirNovac = fullfile(dirNovac,fileNameStrNovac); % Determines full directory to access data 
    [asp] = plotTimeseries(inDirNovac, vName, formatInNovac, fileTypeGvp, dirGvp, formatInGVP, colSO2, minWidth, colVEI, x0, y0, wX, hY, outDirFig, outDirMat, asp);
 
    subplot(2,1,2)
    fileName = volcano(volIndx); vName = char(fileName); vn = string(vName); % Determines filename according to volcano
    fileNameStrMat1 = [vName ext2bs fileTypeMat]; exOut1 = exist(fileNameStrMat1, 'file'); 
    fileNameStrMat2 = [vName ext2nbs fileTypeMat]; exOut2 = exist(fileNameStrMat2, 'file');
    if exOut1 == 2 && exOut2 == 2
        load(fileNameStrMat1, 'dataSet', 'fMasterBs', 'pthMasterBs', 'msdMaster'); % Bootstrapped parameters
        load(fileNameStrMat2, 'pxxMasterNbs', 'fMasterNbs', 'pthMasterNbs'); % Non-bootstrapped parameters
        plotLsAllPoverlay(vn, fMasterBs, dataSet, fMasterNbs, pthMasterNbs, pthMasterBs, msdMaster, pxxMasterNbs, colSel, deltaT, Nf, colSD, lim, x0, y0, wX, hY, outDirFig, Pfa);
    else
    end
end
    
% END OF SCRIPT