function [Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,Pop_casual,Pop_casual_history,Pop_desire,num_born,IdCounter]=Birth(Pop_Id,Pop_Age,Pop_Steady,Pop_steady_history,Pop_casual,Pop_casual_history,Pop_desire,CasualProp,lambda_d,IdCounter)
% give birth to the new lambda_y individuals
% set their age to 15 years, set their steady and casual list of partners
% to be empty, set the respective history entries to -1, seed desire for
% casual partner acquisition
num_born=poissrnd(lambda_d);
if num_born>0
    ind=find(Pop_Id==0,num_born);
    Pop_Age(ind)=15;
    new_ids=IdCounter:1:(IdCounter+num_born-1);
    Pop_Id(ind)=new_ids;
    Pop_Steady(1:2,ind)=0;
    Pop_casual(:,ind)=0;
    Pop_steady_history(:,ind)=-1;
    Pop_casual_history(:,ind)=-1;
    IdCounter=IdCounter+num_born;
    x=0:1:50; % prefered number of casual partners within 12 months

    cfX=cumsum(CasualProp(1,:));
    for counter=1:1:num_born
        r1=rand;
        ind0=find(cfX>r1,1);
        Pop_desire(1,ind(counter))=x(ind0);%desired number of partners within the last 12 months
        Pop_desire(2,ind(counter))=0;
    end
end
end
