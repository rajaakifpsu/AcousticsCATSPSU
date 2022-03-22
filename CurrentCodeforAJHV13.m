%v13: changed the whole thing to a function taking filepaths and savepath
%as parameters

%v12: change font size
%corrected filepaths
%functions for displaying plots changed, but not any calculations.

%Purpose: to analyze acoustics data, spit out a ton of plots

%inputs:
%       filename is string, complete file name and location for the
%       acoustics data, a .h5 file. example: filename = "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\NASA Tmotor\CCWCases\Steady\tdhcs7\acs_data.h5";
%       perFileName is that but for the performance data, a .txt file
%       noiseName is the string name of the .mat file where data for
%       background noise with psu, background noise without power supply
%       unit, and motor noise are all stored. (They are calculated in another code so as not to waste time repeating calculations)
%       isLong is a boolean for if the test is long (120s). Needed so the
%       code can select 30s of data. If the code is short, it is already
%       30s of data so it doesnt need to select.
%       savePath is the complete name and location of the folder where
%       graphs should be saved to.


function CurrentCodeforAJHV13(filename,perfFileName,noiseName,motorNoiseFilePath, isLong, savePath)

%% Loading in data for performance and acoustics

tic

numofmics = 12;
HPend=100; %where the highpass filter will go to (Hz)

file = hdf5info(filename); %the acoustics data
performanceMatrix=readmatrix(perfFileName); %the performance data

load(noiseName) %the background noise (with and without psu). They are 
%calculated elsewhere in another code, CurrentCodeforBackgroundNoise.m

PSDMotorNoiseRed=CurrentCodeforMotorNoise(HPend, numofmics, motorNoiseFilePath); %PSD of motor noise

%Reading all the Signals from the h5 file
%to work with the processing techniques
fs = hdf5read(file.GroupHierarchy.Datasets(4)); %sampling rate in Hz. Which is the 4th field
acousticsdata_uncorr = hdf5read(file.GroupHierarchy.Datasets(1)); % Acoustics value in voltage.

%% Downsampling the signals

acousticsdata_uncorr = acousticsdata_uncorr(1:2:length(acousticsdata_uncorr(:,1)),:); %downsampling the signals for computational efficiency
fs=fs/2; %downsampling the signal for computational efficiency

%% Converting from voltage to Pa with microphone sensitivities

micsensitivities = hdf5read(file.GroupHierarchy.Datasets(5)); % Microphone sensitivities mv/Pa
acousticsdata = (1000./(micsensitivities')).*acousticsdata_uncorr; % Calibrated Acoustic values%initialize calculations
noofsamples = length(acousticsdata); %length of signal

%% Correct the distance from rotor to microphones to be the distance NASA wants.

correctDist=1.8956; %in m (the distance NASA wants)
micDist=[65.19, 62.97, 61.34, 60.34, 60.00, 60.34, 61.34, 62.97, 65.19, 67.93, 71.14, 74.75];%in inches. distance we have to each mic in test.
micDist=convlength(micDist, 'in','m'); %converts micDist to m

for i=1:numofmics
    acousticsdata(:,i) = acousticsdata(:,i) * (micDist(i)/correctDist); %its just signal multiplied by (dist we have / dist we want)
end

%% Defining other little variables

N = noofsamples;
dt = 1/fs;    %time increment between each data point
time_array = [0:N-1]*dt;  %get array of times at each point
T = max(time_array); %max time or period
df = 1/T; 

%% take all the pressure time data here

for i = 1:numofmics
    acoustics_total(:,i) = acousticsdata(:,i); 
    pressure_total(:,i) = acoustics_total(:,i); 
end 

%%  taking a section of time in the .h5 file to define pressure_reduced

starttime = 90; %start time in seconds
endtime = 120; %end time in seconds
T_reduced = endtime - starttime;
N_reduced = T_reduced*fs;

if isLong %if the data is long, select 30s of it
    for i = 1:numofmics %for every microphone
        acoustics_reduced(:,i) = acousticsdata((starttime*fs:endtime*fs-1),i); %sound pressures
        time_array_reduced(:,i) = time_array(starttime*fs:endtime*fs-1); %getting the reduced time for each mic
        pressure_reduced(:,i) = acoustics_reduced(:,i);
    end
else %if the testis already 30s long
    acoustics_reduced=acousticsdata;
    time_array_reduced=time_array;
    pressure_reduced=acoustics_reduced;
end

%% Correct the DC Offset

for i = 1:numofmics
   
   meanDCoffset(:,i) = mean(pressure_total(:,i),1);
   pressure_total(:,i) = pressure_total(:,i) - meanDCoffset(1,i);
   pressure_reduced(:,i) = pressure_reduced(:,i) - meanDCoffset(1,i);
    
end

%% calculate the OASPL RMS

 dbref = 2e-5; % reference pressure in Pa or 2e-5;
 P_ref = dbref;
 G_ref = dbref^2; 

for k = 1:numofmics
    rmspress_total(k,:) = rms(pressure_total(1:N,k));
    OASPL_total_RMS(k,:) = 20*log10(rmspress_total(k,:)/P_ref); %OASPL at each Mic
    OASPL_reduced_RMS(k,:) = 20*log10(rms(pressure_reduced(1:N_reduced,k))/P_ref); %OASPL reduced at each Mic
    deltatotalreduced = OASPL_total_RMS - OASPL_reduced_RMS;
end 


%% Angles of the Microphone Array

angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert angles to radians

horzdistance_centerMic = 60*0.0254; %inches to meters
vertdistance_centerMic = 78*0.0254; %inches to meters

%% define RPM
performanceSampRate=20000;
RPM = getRPMSteady(performanceMatrix, performanceSampRate);
numblades = 2;
%% P welch method
 bandwidth = 20; % width of the record in frequency (Hz)
 binwidth = 2^nextpow2(fs/bandwidth); %not necesary for math, but reduces time needed to run code.

 %doing pwelch on the full set of data
%%%[pxxTot,fpwelchTot, PSDTot,PowerSpectrumTot, oasplTot, pxxBBTot, PSDBBTot, PowerSpectrumBBTot, oasplBBTot, pxxVKTot, PSDVKTot, PowerSpectrumVKTot, oasplVKTot, filteredTot, residualBBTot] = pwelchFullAndBB(pressure_total, binwidth,fs,P_ref,numofmics, RPM, numblades);

%now on the reduced set of data
[pxxRed,fpwelchRed, PSDRed,PowerSpectrumRed, oasplRed, spl100Red, spl500Red, spl1000Red, spl10kRed, pxxBBRed, PSDBBRed, PSDBBunsmoothedRed, PowerSpectrumBBRed, oasplBBRed, spl100BBRed, spl500BBRed, spl1000BBRed, spl10kBBRed, pxxVKRed, PSDVKRed, PowerSpectrumVKRed, oasplVKRed, spl100VKRed, spl500VKRed, spl1000VKRed, spl10kVKRed, pxxHPORed, PSDHPORed, PowerSpectrumHPORed, oasplHPORed, spl100HPORed, spl500HPORed, spl1000HPORed, spl10kHPORed, filteredRed, residualBBRed, highpassOnlyRed] = pwelchFullAndBB(pressure_reduced, binwidth,fs,P_ref,numofmics, RPM, HPend, numblades);
                                                   
%% Plots
%Graphs power spectrum aka narrowband spectrum vs frequency (all 12 mics)
% % % graphTitlesTot=["Total PSD for the whole range of points","Total PSD for the whole range of points","Total PSD for microphone "];
% % % createPowerSpectralDensityPlots(fpwelchTot, PSDTot, numofmics, graphTitlesTot);
% % % 
% % % graphTitlesBB=["Broadband PSD for the whole range of points","Broadband PSD for the whole range of points","Broadband PSD for microphone "];
% % % createPowerSpectralDensityPlots(fpwelchBBTot, PSDBBTot, numofmics, graphTitlesBB);

set(0,'defaultAxesFontSize',16); %sets default font size to 16

Setcol=["--b","b","g","--k","k","--r","r","m"]; %for formatting the line types in PSD plots of different data sets. 
%Order: Unfiltered, Highpassonly, Vold-Kalman,Broadband (with Peaks),Broadband (without Peaks),Background Noise without PSU,Background Noise with PSU,Motor Noise

MicCol=[135, 0, 0; 135, 77, 0; 135, 126, 0; 119, 135, 0; 74, 135, 0; 0, 135, 104; 0, 126, 135; 0, 74, 135; 0, 14, 135; 59, 0, 135; 112, 0, 135; 135, 0, 97;]/255;
%for formatting the line types in PSD plots of different mic numbers. RGB triplets used. 
%Order: mics 1-12 in order


graphTitlesRed=["Unfiltered PSD","Unfiltered PSD","Unfiltered PSD for microphone "];
createPowerSpectralDensityPlots(fpwelchRed, PSDRed, numofmics, numblades, RPM(2), graphTitlesRed, 0, MicCol, savePath);
saveDeletePlots(savePath);

graphTitlesBBRed=["Broadband PSD","Broadband PSD","Broadband PSD for microphone "];
createPowerSpectralDensityPlots(fpwelchRed, PSDBBRed, numofmics, numblades, RPM(2), graphTitlesBBRed, HPend, MicCol, savePath);
saveDeletePlots(savePath);

%Directivity
angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
angle = angle.*(pi/180); %convert into radians (only for plotting)
%%%createDirectivityPlot(oasplTot,oasplBBTot,numofmics,'Full-Data Directivity Plot')
directivityTitles=["Directivity Plot of OASPL [0.1 20]kHz","Directivity Plot of SPL in range [0.1, 20] kHz","Directivity Plot of SPL in range [0.5, 20] kHz","Directivity Plot of SPL in range [1, 20] kHz","Directivity Plot of SPL in range [10, 20] kHz","Directivity Plot", "Broadband Directivity Plot"];
createDirectivityPlot(oasplRed,oasplBBRed,oasplHPORed,spl100Red,spl100BBRed,spl100HPORed,spl500Red,spl500BBRed,spl500HPORed,spl1000Red,spl1000BBRed,spl1000HPORed,spl10kRed, spl10kBBRed, spl10kHPORed, pxxBBRed, fpwelchRed, P_ref, fpwelchRed(2)-fpwelchRed(1) ,numofmics,directivityTitles)

saveDeletePlots(savePath);

%Directivity of BPF and Harmonics therof
% % % SPLForBPFHarmonicsTot=createDirectivityPlotForHarmonics(pxxTot,oasplTot,fpwelchTot,RPM(2),numblades,numofmics,P_ref,'BPF Harmonics SPL Directivity Plot (Total Data)');
SPLForBPFHarmonicsRed=createDirectivityPlotForHarmonics(pxxRed,spl100Red,fpwelchRed,RPM(2),numblades,numofmics,P_ref,'BPF Harmonics SPL Directivity Plot (Reduced Data)');
saveDeletePlots(savePath);

%Graphs Greenwood asked for to see if VK works
%createCheckGraphs(filteredTot, residualBBTot, pxxTot, PSDTot, PSDBBTot, P_ref, fs, binwidth, numblades, numofmics, RPM)
PSDatBPFHarmonics=createCheckGraphs(fpwelchRed, pressure_reduced, filteredRed, residualBBRed, pxxRed, PSDRed, PSDBBRed, PSDBBunsmoothedRed, PSDVKRed, PSDHPORed, PSDBGnoPSURed, PSDBGyesPSURed, PSDMotorNoiseRed, P_ref, fs, binwidth, numblades, numofmics, RPM, HPend, Setcol, MicCol, savePath);

saveDeletePlots(savePath);

%Graphs for NASA
% % % createNasaDeliverables(fpwelchTot, pxxTot,SPLForBPFHarmonicsTot, P_ref, angle, numofmics);
createNasaDeliverables(SPLForBPFHarmonicsRed, angle, numofmics, spl500Red, oasplRed, spl100Red);
saveDeletePlots(savePath);


toc






%% Extra stuff- junk code but dont delete

% PERFORMANCE STUFF - will require more tests to be worth something.
%====================================================================
% sampleRate=20000; %samplerate of the performance measurments
% tempTdhcs=13.4; %temperture (C) for the tdhcs test 7. Manually copied from TestMatrixRaja_SepEntry_II
% pressTdhcs=976;%same but for pressures (mbar)
% testDataTdhcs=performanceMatrix;
% testDataTdhcsCorr=densityCorrection(testDataTdhcs,tempTdhcs,pressTdhcs);
% 
% [TDHCSavFMRPM,TDHCSavError,TDHCSstanDevs,TDHCSSEM,~,~,~,~]=getDataFromWhole(testDataTdhcs, sampleRate,false);
% [TDHCSavFMRPMCorr,TDHCSavErrorCorr,TDHCSstanDevsCorr,TDHCSSEMCorr,~,~,~,~]=getDataFromWhole(testDataTdhcsCorr, sampleRate,false);
% 
% createFnMProfiles1(TDHCSavFMRPM,TDHCSSEM,TDHCSavFMRPMCorr,TDHCSSEMCorr,100:101)
% createFnMProfiles2n3(TDHCSavFMRPM,TDHCSSEM,TDHCSavFMRPMCorr,TDHCSSEMCorr,102:103,[0,15;-.25,0.05])










%%%A BUNCH OF STUFF I CAN PROBABLY DELETE

% % % 
% % % for k = 1:numofmics
% % % 
% % %     fq = 0:bandwidth:fs/2;
% % %     %[pxx,f] = pwelch(pressure_total(:,k) ,W,nooverlap,[],fs);
% % %     %[pxx(:,k),f(:,k)] = pwelch(pressure_total(:,k),W,nooverlap,[],fs);
% % %     %W is a vector
% % %     [pxx(:,k),fq(:,k)] = pwelch(pressure_total(:,k),80000,512,[],fs); %80000 is interget (where issue occurs)
% % %     %[pxx(:,k),f(:,k)] = pwelch(pressure_reduced(:,k),W,nooverlap,[],fs);
% % %     
% % %     
% % %     %[pxx,f] = pwelch(x,W,nooverlap,f,fs) returns the two-sided Welch PSD estimates at the frequencies specified in the vector, f. The vector f must contain at least two elements, because otherwise the function interprets it as nfft. The frequencies in f are in cycles per unit time. The sample rate, fs, is the number of samples per unit time. If the unit of time is seconds, then f is in cycles/sec (Hz).
% % %     
% % %     %from miccal
% % % %     Spl_denorm = pxx / bandwidth; 
% % % %     energy_tot = sum(Spl_denorm);
% % % %     oaspl = 10*log10((energy_tot)/Pref.^2);
% % % 
% % %     
% % %     %SPL(:,i) = 10*log10(Gxx(:,i)./Gref);
% % %     %SPLp(:,k) = 10*log10(pxx(:,k)./G_ref); %SPL at each Mic
% % %     SPLp(:,k) = 10*log10(pxx(:,k)./G_ref/bandwidth); %SPL at each Mic
% % % 
% % %     %valp(:,k) = (1/T)*(sum((pxx(:,k)))); %energy in the system
% % %     %OASPLp(:,k) = 10*log10(valp(:,k)./G_ref); %OASPL at each Mic  
% % %    
% % %     
% % %     Spl_denorm(:,k) = pxx(:,k) / bandwidth; 
% % %     energy_tot(:,k) = sum(Spl_denorm(:,k));
% % %     oaspl(:,k) = 10*log10((energy_tot(:,k))/P_ref.^2);
% % %     %RMS_Gxx(:,k) = sqrt(sum((Gxx_use(:,k)))*df); %root mean squared energy. or root of energy.
% % %     %OASPL_exp_gxxmethod(:,k) = 20*log10(RMS_Gxx(:,k)./P_ref);
% % %     
% % %     RMS_pxx(k) = sqrt(sum(pxx(:,k))*df); %root mean squared energy. or root of energy.
% % %     OASPL_exp_pxxmethod(k) = 10*log10(RMS_pxx(k)./G_ref);
% % % 
% % % 
% % %     figure()
% % %     hh = plot(f(1:nooverlap-1,k), pxx(1:nooverlap-1,k), '-b'); %skip point 1
% % %     set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% % %     set(gca, 'fontsize', 13) % gca = "get current axis"
% % %     grid on
% % %     xlim([0 20000]) % x limit to half of samplerate
% % %     title({'Gxx for mic ',num2str(k)})
% % %     xlabel('Frequency [Hz]','fontsize', 13)
% % %     ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)
% % % 
% % %     figure()
% % %     hh = plot(f(61:nooverlap-1,k), SPLp(61:nooverlap-1,k), '-b'); %skip point 1
% % %     set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% % %     set(gca, 'fontsize', 13) % gca = "get current axis"
% % %     grid on
% % %     xlim([0 20000]) % x limit to half of samplerate
% % %     title({'SPL for mic ',num2str(k)})
% % %     xlabel('Frequency [Hz]','fontsize', 13)
% % %     ylabel('SPL [dB]', 'fontsize', 13)
% % % 
% % % end
% % %     
% % %     angle = [23.01, 17.67, 11.99, 6.06, 0.00, -6.06, -11.99, -17.67, -23.01, -27.96, -32.50, -36.6];
% % %     angle = angle.*(pi/180);
% % %     figure()
% % %     rlim = [30 45]; %dBs
% % %     tlim = [-36.62 23.01];  %degrees
% % %     polarplot(angle(1:11),OASPLp(1:11),'LineWidth',1.3)
% % %     ax=gca;
% % %     ax.FontSize=14;
% % %     grid on
% % %     grid minor
% % %     ax.RLim=rlim;
% % %     ax.ThetaLim=tlim;
% % %     legend({'w/Filter','Original'})
% % %     title(['OASPL Directivity for ',num2str(RPM),...
% % %           'rpm [Hover]'])
% % %     % saveas(gca,['BLtd Plot for Sep ',num2str(RPM),'rpm[Edgewise].png'])
% % %     % saveas(gca,['Tonal Directivity for Sep ',num2str(RPM),'rpm[Edgewise].png'])
% % % 
% % %     
%% For each Mic to check the value and learn without using any filter or overlaps
% 
% freq_axis = ((1:N)*df)';
%  
% lsp = fft(acoustics_total)*dt; %Find the linear spectrum, lsp 
% Sxx = abs(lsp).^2./T; %get dual sided spectrum,[Pa^2/Hz]. Here Sxx size is numofsamples
% 
%set requirements
% %at N/2| Gxx(N/2) = Sxx(N/2)
% %at 0| Gxx(1) = Sxx(1)
% %others| Gxx = 2*Sxx from 1 to N/2 - 1 in freq domain
% 
% %for mic 1
% mic = 1;
% Sxx_mic = Sxx(:,mic);
% Gxx = 2*Sxx_mic;  %(2/T)*abs(lsp).^2
% Gxx(1) = Sxx_mic(1);
% Gxx(N/2+1) = Sxx_mic(N/2+1); 
% %Gxx has Noofsamples size
% %Sxx has Noofsamples size
% 
% % Compare Mean of Time Series with the convention
% % conventional must = to Ms time series
% 
% % time domain version. Mean-squared value of the original time series.
% MS_timeSeries = sum(acoustics_total(:,1).^2)./N;
% 
% % Sxx*df
% MS_Xm = sum(Sxx_mic*df); % integral of double sided. conventional
% 
% % Gxx*df but Gxx is to half of number of samples
% MS_Gxx = sum(Gxx(1:(N/2+1))*df); %integral of single sided. conventional
% 
% % if the MS time series =/= Sxx*df to all number of samples, Sxx is wrong
% % if the MS time series =/= Gxx*df to half number of samples, Gxx is wrong
% 
% % plot the Gxx to half of the number of samples
% figure()
% hh = plot(freq_axis(2:(N/2)+1), Gxx(2:(N/2)+1), '-b'); %skip point 1
% set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% set(gca, 'fontsize', 13) % gca = "get current axis"
% grid on
% xlim([0 40000]) % x limit to half of samplerate
% title({'Gxx for mic ',num2str(mic)})
% xlabel('Frequency [Hz]','fontsize', 13)
% ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)
% 
% %get log plot
% figure()
% hh = semilogy(freq_axis(2:(N/2)), Gxx(2:(N/2)+1), '-b'); %skip point 1
% set(hh, 'linewidth', 1) % hh is "handle" for "plot"
% set(gca, 'fontsize', 13) % gca = "get current axis"
% xlim([0 40000]) % x limit to half of samplerate
% grid on
% title({'Log plot of Gxx for mic ',num2str(mic)})
% xlabel('Frequency [Hz]','fontsize', 13)
% ylabel('Gxx [Pa^{2}/Hz]', 'fontsize', 13)

%%







end