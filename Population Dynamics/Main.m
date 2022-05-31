% This script generates population dynamics in the population of men who
% have sex with men. The following processes are included: birth,
% background death, aging out of the population of interest, creation and
% dissolution of steady and casual relationships.
% A single trajectory is generated and is stored in .csv format. The following
% are stored as time series population size, size of age compartments, duration of steady and 
% casual relationships, number of individuals available for a steady 
% relationship, number of individuals in 0,1 or 2 steady relationships,
% number of individuals in 0, 1, 2 , 3, 4, 5 and more than 5 casual and 
% steady relationships in the last XX months. End of the run snapshot:
% Exact age distribution of the population, distribution of the current
% number of casual partners and steady partners (at the end and within the
% last XX months).
%%
% preparation of state space
clear variables;
close all;
clc;
format long;
myname=mfilename;
error_tol=1e-6;
edges1year=15:1:75;
edges10year=15:10:75;
num_age=numel(edges10year)-1;
year=365;% duration of a year

%%
% reading of the distributions and parameters
FolderStr='Data/';
AgeFileStr='AgeDistr.csv'; % array of .csv setting distribution between the ages of
% of 15 and 75 by increments of 10, should add up to 1
DeathRateFileStr='ProbDeathAge.csv';
PerCurSteadyStr='NumberCurrentSteady';
CasualPropAgeStr='CasualDesDistrAge.csv';% distribution of propensity to acquire casual relationships
ParFileStr='Parameters.txt';% file containing lose parameters

AgeDistr=readmatrix([FolderStr,AgeFileStr]);
DeathRates=readmatrix([FolderStr,DeathRateFileStr]);
PerCurSteady=readmatrix([FolderStr,PerCurSteadyStr]);
CasualProp=readmatrix([FolderStr,CasualPropAgeStr]);

% read the parameter file
fid=fopen([FolderStr,ParFileStr],'rt');%open files in text mode
while true
  thisline = fgetl(fid);
  if ~ischar(thisline)
      break;
  end  %end of file
  % turn string into command  
  eval(thisline);
end
fclose(fid);

% validating the input
if sum(AgeDistr<0)>0
    error([myname,': age distribution should have non-negative values']);
end

if abs(sum(AgeDistr)-1)>error_tol
    error([myname,': age distribution should add up to one']);
end

if ~exist("T",'var')
    T=1;
    disp([myname,': the total duration of simulation was not defined, was set to ',num2str(T)]);
end

if ~exist("lambda",'var')
    error([myname,': birth rate was not supplied']);
end

lambda_d=lambda/year;
IdCounter=N+1;
%%
% Allocating the arrays and setting the population to its initial state
% ids
Pop_Id=zeros(1,2*N);
%age
Pop_Age=zeros(1,2*N);
%steady relationships
Pop_Steady=zeros(2,2*N);
Pop_steady_history=zeros(Max_parts,2*N)-1;% history of acquisition of steady partners within the last 
Steady_dur=[];

% casual relationships
CasualRels=[];
%duration of relationships: all of the values will be preserved to obtain
%the distribution
Casual_dur=[];
% container for the current number of partners
Pop_casual=zeros(Max_parts,2*N);% the first row is for ids
%container for the history in the last XX months
Pop_casual_history=zeros(2*Max_parts,2*N)-1;% the first row is for ids

Pop_Id(1,1:1:N)=1:1:N;
% assign age
%allocate the container to track time series of age compartments
yAge=zeros(num_age,T*year+1);
% age brackets: 15-25, 25-35, 35-45, 45-55, 55-65, 65-75
CSAge=cumsum(AgeDistr);
% distribute the propensity to acquire casual partnership
Pop_desire=zeros(2,2*N);
x=0:1:50; % prefered number of casual partners within 12 months
Pop_summar_desire_age=zeros(num_age,numel(x));
for counter=1:1:N
    r1=rand;
    ind=find(CSAge>r1,1);
    yAge(ind,1)=yAge(ind,1)+1;
    a1=15+(ind-1)*10;
    Pop_Age(1,counter)=unidrnd(10)+a1-1;
    % based on the age of the individual determine their propensity to
    % acquire casual partners
    cfX=cumsum(CasualProp(ind,:));
    r2=rand;
    ind2=find(cfX>r2,1);
    Pop_desire(1,counter)=x(ind2);%desired number of partners within the last 12 months
    Pop_summar_desire_age(ind,x(ind2)+1)=Pop_summar_desire_age(ind,x(ind2)+1)+1;
end

for counter=1:1:num_age
    Pop_summar_desire_age(counter,:)=Pop_summar_desire_age(counter,:)/sum(Pop_summar_desire_age(counter,:));
end

% output check for the of desire to form partnerships vs the data
figc=1;
figure(figc);
subplot(1,2,1);
bar(CasualProp,'stacked');
ylim([0,0.3]);
xlabel('Age','interpreter','latex');
ylabel('Probability','interpreter','latex');
title('Source distribution','interpreter','latex');

%rewrite ticks
xticks(1:1:num_age);
xlabels={};
counter0=1;
for counter=1:1:num_age
    str=[num2str(15+(counter-1)*10),'-',num2str(15+counter*10)];
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);
str={};
for counter=0:1:50
    str{counter+1}=num2str(counter);
end
legend(str,'Location','southoutside');
set(gca,'FontSize',25);

subplot(1,2,2);

bar(Pop_summar_desire_age,'stacked');
ylim([0,0.3]);
%rewrite ticks
xticks(1:1:num_age);
xlabels={};
counter0=1;
for counter=1:1:num_age
    str=[num2str(15+(counter-1)*10),'-',num2str(15+counter*10)];
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);

xlabel('Number of partners','interpreter','latex');
ylabel('Proportion');
title('Resulting distribution','interpreter','latex');
set(gca,'FontSize',25);
figc=figc+1;

%generate steady relationships
%allocate the container to track time series of age compartments
ySteady=zeros(3,T*year+1);
% calculate the number of possible new pairs
ind1=find(Pop_Steady(1,1:1:N)==0);
ind2=find(Pop_Steady(2,1:1:N)==0);
ind=sort(unique([ind1,ind2]));
pot_pairs=floor((numel(ind)-1)*numel(ind)/2);
Rels_steady=[];
num_new_pairs=0;
for counter=1:1:pot_pairs
    r1=rand;
    if r1<sigma
        s_fl=0;
        num_iter=0;
        while ~s_fl & num_iter<N_Iter
            [s_fl,Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady]=make_steady_pairV2(Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady,psi,Pop_casual,0);
            num_iter=num_iter+1;
        end
        if s_fl>0
            num_new_pairs=num_new_pairs+1;
        end
    end
end

disp(['Out of ',num2str(pot_pairs),' potential new pairs, ',num2str(num_new_pairs),' were created']);

% summarize the population statistics for the output
n_single=numel(find(Pop_Steady(1,1:1:N)==0));
n_atleastone=N-n_single;
n_two=numel(find(Pop_Steady(2,1:1:N)~=0));
n_one=n_atleastone-n_two;
ySteady(1,1)=n_single;
ySteady(2,1)=n_one;
ySteady(3,1)=n_two;

ySteady12=ySteady;

%generate casual relationships
%allocate the container to track time series of age compartments
yCasual=zeros(7,T*year+1);%0, 1 , 2, 3, 4, 5, and more than 5
% calculate the number of possible new pairs

Rels_casual=[];
num_new_pairs=0;

for counter=1:1:(N*(N-1)/2)
    r1=rand;
    if r1<sigma2
        s_fl=0;
        num_iter=0;
        while ~s_fl & num_iter<N_Iter
            [s_fl,Pop_casual,Pop_casual_history,Pop_desire,Rels_casual]=create_casual_pairV2(0,Pop_Id,Pop_casual, Pop_casual_history,Rels_casual,Pop_desire,Pop_Steady);
            num_iter=num_iter+1;
        end
        if s_fl>0
            num_new_pairs=num_new_pairs+1;
        end
    end
end

% summarize the population statistics for the output
yCasual(1,1)=numel(find(Pop_desire(2,1:N)==0));
yCasual(2,1)=numel(find(Pop_desire(2,1:N)==1));
yCasual(3,1)=numel(find(Pop_desire(2,1:N)==2));
yCasual(4,1)=numel(find(Pop_desire(2,1:N)==3));
yCasual(5,1)=numel(find(Pop_desire(2,1:N)==4));
yCasual(6,1)=numel(find(Pop_desire(2,1:N)==5));
yCasual(7,1)=N-sum(yCasual(:,1));
yCasual12=yCasual;

% Demographics initial settings

figure(figc);
bar(AgeDistr,'BarWidth',1);
xlabel('Age, years','Interpreter','latex');
ylabel('Proportion','Interpreter','latex');
title({'Source age distribution'},'Interpreter','latex');
xticks(1:1:num_age)
xticklabels({'15-25','25-35','35-45','45-55','55-65','65-75'});
set(gca,'FontSize',25);
figc=figc+1;

figure(figc);
histogram(Pop_Age(1:1:N),edges1year,'Normalization','probability');% histogram with step of 1 year
xlabel('Age, years','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
title({'Age distribution at';'the start of the simulaiton'},'Interpreter','latex');
set(gca,'FontSize',25);
figc=figc+1;

figure(figc);
h1=histogram(Pop_Age(1:1:N),edges10year,'Normalization','probability');% histogram with step of 1 year
xlabel('Age, years','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
title({'Age distribution at';'the start of the simulaiton'},'Interpreter','latex');
xticks(20:10:70);
xticklabels({'15-25','25-35','35-45','45-55','55-65','65-75'});
set(gca,'FontSize',25);
figc=figc+1;

% Sexual networks initial settings
% steady relationships
figure(figc);
bar(PerCurSteady);
xlabel('Number of current steady partners','Interpreter','latex');
ylabel('Proportion','Interpreter','latex');
title({'Source current steady partner distribution'},'Interpreter','latex');
xticks(1:1:3)
xticklabels({'Single','One partner','More than one partner'});
set(gca,'FontSize',25);
figc=figc+1;

% % casual relationships
% figure(figc);
% plot(x,CasualProp,'ro','MarkerSize',6);hold on;
% xlim([0,max(x)]);
% ylim([0,0.2]);
% xlabel({'Number of casual';'partners within 12 months'},Interpreter='latex');
% ylabel('Probability','interpreter','latex');
% title('The source distribution','Interpreter','latex');
% set(gca,'FontSize',25);
% figc=figc+1;

%%
% find ids of these who are available
indexAlive=find(Pop_Id>0);
ind_avail1=find(Pop_Steady(1,indexAlive)==0);
ind_avail2=find(Pop_Steady(2,indexAlive)==0);
ind_avail=sort(unique([ind_avail1,ind_avail2]));

% Main simulation
for tcounter=1:1:(365*T)
    % ageing module - shift the age
    if mod(tcounter,365)==0
        Pop_Age=Age(Pop_Age,Pop_desire,CasualProp);
    end
    % death module 
    [Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,Rels_steady,Steady_dur,Pop_casual,Pop_casual_history,Rels_casual,Casual_dur,Pop_desire]=Death(Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,DeathRates,Rels_steady,Steady_dur,Pop_casual,Pop_casual_history,Rels_casual,Casual_dur,Pop_desire,tcounter); 
    % steady relationship creation
    
    indexAlive=find(Pop_Id>0);
    ind_avail1=find(Pop_Steady(1,indexAlive)==0);
    ind_avail2=find(Pop_Steady(2,indexAlive)==0);
    candidates_ids=sort(unique([ind_avail1,ind_avail2]));

    %create a list of possible new pairs
    Pot_pairs=CreatePairList(candidates_ids);
    if size(Rels_steady,1)>0
        Rels_steady_t=sortrows(Rels_steady);
        ind_ineligible=[];
        for counter=1:1:size(Rels_steady_t,1)
            [flag,index]=ismember(Rels_steady_t(counter,1:2),Pot_pairs,'rows'); %remove from potential candidates these pairs that already exist
            if flag
                Pot_pairs(index,:)=[];
            end
        end
    end
    
    % take a note of the size of the Rels array to know which relationships
    % were just added as not to accidentally break them up
    num_new_pairs=0;
    for counter=1:1:size(Pot_pairs,1)
        r1=rand;
        if r1<sigma_d % create the pair
            [s_fl,Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady]=make_steady_pairV2(Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady,psi,Pop_casual,tcounter);
            num_new_pairs=num_new_pairs+1;
        end
    end

    % steady relationship disolution
    exist_pairs_num=size(Rels_steady,1)-num_new_pairs;
    pairs_list=[];
    for counter=1:1:exist_pairs_num
        r1=rand;
        if r1<rho_d % dissolve the pair
            pairs_list=[pairs_list; counter];
        end    
    end 
    if numel(pairs_list)>0
        % split all the steady couples that should be split
        [Pop_Steady,Rels_steady,Steady_dur]=split_steady_pair_list(Pop_Id,Pop_Steady,Rels_steady,pairs_list,tcounter,Steady_dur);
    end
    %% casual relationships
    % creation
    
    candidates_ids=find(Pop_Id>0);
    
    %create a list of possible new pairs
    Pot_pairs=CreatePairList(candidates_ids);

    %check which of these created already exist as steady partnerships
    if size(Rels_steady,1)>0
        Rels_steady_t=sortrows(Rels_steady);
        ind_ineligible=[];
        for counter=1:1:size(Rels_steady_t,1)
            [flag,index]=ismember(Rels_steady_t(counter,1:2),Pot_pairs,'rows'); %remove from potential candidates these pairs that already exist
            if flag
                Pot_pairs(index,:)=[];
            end
        end
    end

    %check which of these created already exist as casual partnerships
    if size(Rels_casual,1)>0
        Rels_casual_t=sortrows(Rels_casual);
        ind_ineligible=[];
        for counter=1:1:size(Rels_casual_t,1)
            [flag,index]=ismember(Rels_casual_t(counter,1:2),Pot_pairs,'rows'); %remove from potential candidates these pairs that already exist
            if flag
                Pot_pairs(index,:)=[];
            end
        end
    end
    
    % take a note of the size of the Rels array to know which relationships
    % were just added as not to accidentally break them up
    num_new_pairs=0;
    for counter=1:1:size(Pot_pairs,1)
        r1=rand;
        if r1<sigma2_d % create the pair
            [s_fl,Pop_casual,Pop_casual_history,Pop_desire,Rels_casual]=create_casual_pairV2(tcounter,Pop_Id,Pop_casual, Pop_casual_history,Rels_casual,Pop_desire,Pop_Steady);
            num_new_pairs=num_new_pairs+1;
        end
    end
    
    % casual pairs relationship disolution
    exist_pairs_num=size(Rels_casual,1)-num_new_pairs;
    pairs_list=[];
    for counter=1:1:exist_pairs_num
        r1=rand;
        if r1<rho2_d % dissolve the pair
            pairs_list=[pairs_list; counter];
        end    
    end 
    if numel(pairs_list)>0
        % split all the casual couples that should be split
        [Pop_casual,Rels_casual,Casual_dur,Pop_desire]=split_casual_pair_list(Pop_Id,Pop_casual,Rels_casual,Casual_dur,pairs_list,Pop_desire,tcounter);
    end

    %% birth module 
    [Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,Pop_casual,Pop_casual_history,Pop_desire,num_born,IdCounter]=Birth(Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,Pop_casual,Pop_casual_history,Pop_desire,CasualProp,lambda_d,IdCounter);
    
    % tabulate the population, later rewrite to do it only at a handful of
    % points
    ind_active=find(Pop_Age>0);
    index = discretize(Pop_Age(ind_active),edges10year);
    for counter=1:1:num_age
        yAge(counter,tcounter+1)=sum(index==counter);
    end

    %collect the time series data for steady relationships
    ind1=find(Pop_Steady(1,ind_active)==0);
    ind2=find(Pop_Steady(2,ind_active)==0);
    ind=intersect(ind1,ind2);
    n_single=numel(ind);
    n_atleastone=numel(ind_active)-n_single;
    
    ind1=find(Pop_Steady(1,ind_active)>0);
    ind2=find(Pop_Steady(2,ind_active)>0);
    ind=sort(unique([ind1,ind2]));
    n_two=numel(find(ind));
    
    n_one=n_atleastone-n_two;
    ySteady(1,tcounter+1)=n_single;
    ySteady(2,tcounter+1)=n_one;
    ySteady(3,tcounter+1)=n_two;

    % clean up number of partners acquired within
    % the last 12 months
    ind_alive=find(Pop_Id>0);

    for ind=ind_alive % iterate through all individuals in the population
        id=Pop_Id(1,ind);
        % steady
        ind_old=find((tcounter-Pop_steady_history(1:end,ind))>hist_dur);
        % remove these in ind_old from People_history
        Pop_steady_history(ind_old,ind)=-1;
        % casual
        ind_old=find((tcounter-Pop_casual_history(1:end,ind))>hist_dur);
        % remove these in ind_old from People_history
        Pop_casual_history(ind_old,ind)=-1;
    end

    % record time series for the current number of partners

    % steady
    yCasual(1,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==0));
    yCasual(2,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==1));
    yCasual(3,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==2));
    yCasual(4,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==3));
    yCasual(5,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==4));
    yCasual(6,tcounter+1)=numel(find(Pop_desire(2,ind_alive)==5));
    yCasual(7,tcounter+1)=numel(ind_alive)-sum(yCasual(:,1));
    
end
%%
% Graphic output and saving of the state space
figure(figc);
plot(1:1:(T*year+1),sum(yAge,1),'LineWidth',3);hold on;
xticks(0:10*365:(T*365+1));
xlabels={};
counter0=1;
for counter=0:10:T
    str=num2str(counter);
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);
xlabel('Time, years','interpreter','latex');
ylabel('Individuals','Interpreter','latex');
title('Total population size','Interpreter','latex');
set(gca,'FontSize',25);
figc=figc+1;

figure(figc);

for counter=1:1:num_age
    plot(1:1:(T*year+1),yAge(counter,:)./sum(yAge,1),'LineWidth',3);hold on;
end
xticks(0:10*365:(T*365+1));
xlabels={};
counter0=1;
for counter=0:10:T
    str=num2str(counter);
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);
xlabel('Time, years','interpreter','latex');
ylabel('Proportion','Interpreter','latex');
legend('15-25','25-35','35-45','45-55','55-65','65-75');
set(gca,'FontSize',25);
figc=figc+1;

%distribution of the age at the end of the run
figure(figc);
ind=find(Pop_Id>0);
histogram(Pop_Age(ind),edges10year,'Normalization','probability');% histogram with step of 1 year
xticks(20:10:70);
xticklabels({'15-25','25-35','35-45','45-55','55-65','65-75'});
set(gca,'FontSize',25);
xlabel('Age, years','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
title({'Age distribution at';'the end of the simulaiton'},'Interpreter','latex');
set(gca,'FontSize',25)
figc=figc+1;

%distribution of the duration of steady relationship
figure(figc);
histogram(Steady_dur/year,'normalization','probability');
xline(mean(Steady_dur/year),'r--','LineWidth',4);
xlabel('Duration of steady relationships, years','interpreter','latex');
ylabel('Probability','interpreter','latex');
set(gca,'FontSize',25)
figc=figc+1;

%time series for the number of steady partners
figure(figc);
for counter=1:1:3
    plot(1:1:(T*year+1),ySteady(counter,:),'LineWidth',3);hold on;
end
xticks(0:10*365:(T*365+1));
xlabels={};
counter0=1;
for counter=0:10:T
    str=num2str(counter);
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);
xlabel('Time, years','interpreter','latex');
ylabel('Number of people','Interpreter','latex');
legend('Single','One steady partner','Two steady partners');
set(gca,'FontSize',25)
figc=figc+1;

% distribution of steady partners at the end of the simulation run

figure(figc);
ind_alive=find(Pop_Id>0);
Num_steady_12=zeros(1,numel(ind_alive));
counter=1;
for ind=ind_alive
    Num_steady_12(counter)=sum(Pop_steady_history(:,ind)>-1);
    counter=counter+1;
end
histogram(Num_steady_12,'normalization','probability');
xlabel(['Number of steady partners within the last',num2str(hist_dur),' months of the simulation'],'interpreter','latex');
ylabel('Probability','interpreter','latex');
set(gca,'FontSize',25);
figc=figc+1;

%% output for casual relationships
%distribution of the duration of casual relationship
figure(figc);
histogram(Casual_dur,'normalization','probability');
xline(mean(Casual_dur),'r--','LineWidth',4);
xlabel('Duration of steady relationships, days','interpreter','latex');
ylabel('Probability','interpreter','latex');
set(gca,'FontSize',25)
figc=figc+1;

%time series for the number of current casual partners
figure(figc);
for counter=1:1:7
    plot(1:1:(T*year+1),yCasual(counter,:),'LineWidth',3);hold on;
end
xticks(0:10*365:(T*365+1));
xlabels={};
counter0=1;
for counter=0:10:T
    str=num2str(counter);
    xlabels{counter0}=str;
    counter0=counter0+1;
end 
xticklabels(xlabels);
xlabel('Time, years','interpreter','latex');
ylabel('Number of people','Interpreter','latex');
legend('0','1','2','3','4','5','>5');
set(gca,'FontSize',25)
figc=figc+1;

figure(figc);
Num_casual_12=zeros(1,numel(ind_alive));
counter=1;
for ind=ind_alive
    Num_casual_12(counter)=sum(Pop_casual_history(:,ind)>-1);
    counter=counter+1;
end
histogram(Num_casual_12,'normalization','probability');
xlabel(['Number of casual partners within the last',num2str(hist_dur),' months of the simulation'],'interpreter','latex');
ylabel('Probability','interpreter','latex');
set(gca,'FontSize',25);
figc=figc+1;