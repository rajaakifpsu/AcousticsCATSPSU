%TO DO: 
%1) For some reason, this code does not work with full filepaths, as
%seen commented out bellow. The full filepath for my folder worked in other
%codes, but the filepath for other folders did not work. This leads me to
%believe its a permission issue or something i dont know. For right now,
%manually copying the files into the folder where this code is run from.
%2) Make it so the motorNoise stuff is calculated allongside
% everything else, not in its own little function.
%3) Represent table created in createCheckGraphs as a figure?
%4) When I increased the font size, this messed up the axes of some of the
%PSD plots (the ones where entire range is scaled by the BPF.) To fix this, maybe
%scale it by 2*BPF instead
%5) Somewhere in this horrid mess of code, 17 unnamed plots are generated with no
%data in them. Find them and remove them.
%6)The PSD for the 2+2 case is too large for the y limits. fix.


%MISC NOTES:
%1) Encountered an odd error where calculateSPL said that rsum was an
%unrecognized function or variable. I can think of no explanation for this.
%It very clearly is not. Before this point, the same function had run dozens of
%times without issue. I reran the code, no issue. I am unsure of what the 
%problem is, then. Fluke?


clear all;
close all ;
clc;

%files correspond to: tdhcs2, tdhcs7, vdhcsbb7

%DO NOT DELETE:
% filenames = ["C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\NASA Tmotor\CCWCases\Steady\tdhcs2\acs_data.h5",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\NASA Tmotor\CCWCases\Steady\tdhcs7\acs_data.h5",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_29_21\vdhcsbb7\acs_data.h5"   
% ]; %acoustics
% 
% perfFileNames=["C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\NASA Tmotor\CCWCases\Steady\tdhcs2\Full.txt",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\NASA Tmotor\CCWCases\Steady\tdhcs7\Full.txt",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_29_21\vdhcsbb7\Full.txt"
% ]; %performance
% 
% motorNoiseFilePath=["C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_15_21\mnb2\acs_data.h5",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_15_21\mnbo7\acs_data.h5",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Runs\11_30_21\mnbbb7\acs_data.h5"
% ]; %motor noise
% 
% savePaths =["C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_12_2022\tdhcs2",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_12_2022\tdhcs7",
%     "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_12_2022\vdhcsbb7"
% ];% the location to save figures to.



filenames = ["acs_data_tdhcs2.h5","acs_data_tdhcs7.h5","acs_data_vdhcsbb7.h5"]; %acoustics

perfFileNames=["tdhcs2.txt","tdhcs7.txt","vdhcsbb7.txt"]; %performance

motorNoiseFilePaths=["acs_data_mnb2.h5","acs_data_mnbo7.h5","acs_data_mnbbb7.h5"]; %motor noise

noiseName='NoisePSDs.mat'; %the background noise (with and without psu).
%They are calculated elsewhere in another code, CurrentCodeforBackgroundNoise.m

isLong= [true true false]; %boolean for if the data is measured over a long time (120s). Used to determine if we need to select 30s of data or not. A false indicates that the entire test is already 30s long.

savePaths =["C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_20_2022\tdhcs2",
    "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_20_2022\tdhcs7",
    "C:\Users\aaron.hafner\OneDrive - The Pennsylvania State University\Share with Acoustics Group\1_Acoustics Code Used\AJHCode\graphs 3_20_2022\vdhcsbb7"
];% the location to save figures to.

for i=length(filenames):-1:1 %analyzes the acoustics for every test specified above. (goes in reverse for debug purposes, doesnt matter)
    CurrentCodeforAJHV13(filenames(i),perfFileNames(i),noiseName,motorNoiseFilePaths(i),isLong(i), savePaths(i));
end