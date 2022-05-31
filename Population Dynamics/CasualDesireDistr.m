% This script generates temporary toy distribution of casual partners
% within the last 12 months

% prepare state space
close all;
clear variables;
clc;
format long;

%destination name string for saving the partner distribution
FolderStr='Data/';
FileStr='CasualDesDistr.csv';
FileStrAge='CasualDesDistrAge.csv';

% potential number of desired partners
x=0:1:50;

% calculate the distribution for the relative rate of partner acquisition
a=5;
b=1/a-0.2;
f=@(x)1./(x+a)-b;

%output the distribution
evalfx=arrayfun(f,x);
evalfx=evalfx/sum(evalfx);
figc=1;
figure(figc);
plot(x,evalfx,'ro','MarkerSize',6);hold on;
xlim([0,max(x)]);
ylim([0,0.2]);

xlabel({'Number of casual';'partners within 12 months'},Interpreter='latex');
ylabel('Probability','interpreter','latex');
title('The source distribution','Interpreter','latex');
set(gca,'FontSize',25)
figc=figc+1;

writematrix(evalfx,[FolderStr,FileStr]);

b_incr_arr=linspace(3,0.3,6);
a=3;
PropensityRates=zeros(6,51);
counter=1;
for b_incr=b_incr_arr
    b=1/30;
    f=@(x)1./(x+a)-b+b_incr;
    
    %output the distribution
    evalfx=arrayfun(f,x);
    figure(figc);
    subplot(1,2,1);
    plot(x,evalfx,'o','MarkerSize',6);hold on;
    xlim([0,max(x)]);
    %ylim([0,0.2]);
    
    xlabel({'Number of casual';'partners within 12 months'},Interpreter='latex');
    ylabel('Non-normalized values','interpreter','latex');
    title('The source distribution','Interpreter','latex');
    set(gca,'FontSize',25);
    
    evalfx=evalfx/sum(evalfx);
    figure(figc);
    subplot(1,2,2);
    plot(x,evalfx,'o','MarkerSize',6);hold on;
    xlim([0,max(x)]);
    %ylim([0,0.2]);
    
    xlabel({'Number of casual';'partners within 12 months'},Interpreter='latex');
    ylabel('Probability','interpreter','latex');
    title('The source distribution','Interpreter','latex');
    set(gca,'FontSize',25);

    PropensityRates(counter,:)=evalfx;
    counter=counter+1;
end
figc=figc+1;
writematrix(PropensityRates,[FolderStr,FileStrAge]);
