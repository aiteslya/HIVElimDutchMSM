% temporary script to wrangle the data. The data is being read from and
% recorded to the folder Data in the root folder.

%%
% Prepare the state space
clear variables;
clc;
close all;
format long;

FolderStr='Data/';

%%
% age distribution
AgeFileSourceStr='MaleNetherlands2020.csv';
AgeFileDestStr='AgeDistr.csv';

ageData=readmatrix([FolderStr,AgeFileSourceStr]);
ageData([1:1:3,16:1:20],:)=[];% remove age brackets 0-15 years and over 75 years
%collapse 5-year age brackets to 10-year age brackets
AgeDataColl=zeros(6,1);
for counter=1:1:6
    AgeDataColl(counter,1)=ageData(2*counter-1,1)+ageData(2*counter,1);
end
%non-dimensionalize the distribution
AgeDataColl=AgeDataColl./sum(AgeDataColl);

writematrix(AgeDataColl,[FolderStr,AgeFileDestStr]);

%% probability of dying between the ages [x; x+5];
ProbDeathFileSourceStr='ProbDyingAgeMaleNetherlands.csv';
ProbDeathFileDestStr='ProbDeathAge.csv';

deathRate=readmatrix([FolderStr,ProbDeathFileSourceStr]);
deathRate(:,1:2)=[];
deathRate(1:2,:)=[];
deathRate=deathRate/(5*365);

deathRateColl=zeros(6,1);
for counter=1:1:6
    deathRateColl(counter,1)=deathRate(2*counter-1,1)+deathRate(2*counter,1);
end

writematrix(deathRateColl,[FolderStr,ProbDeathFileDestStr]);

%% number of current steady partners
NumSteadyPartSourceStr='NumberCurrentSteadyEmis2010.csv';
NumSteadyPartDestStr='NumberCurrentSteady.csv';
NumSteadyParts=readmatrix([FolderStr,NumSteadyPartSourceStr]);
%trim
NumSteadyParts(1,:)=[];
NumSteadyParts(:,1)=[];
writematrix(NumSteadyParts,[FolderStr,NumSteadyPartDestStr]);