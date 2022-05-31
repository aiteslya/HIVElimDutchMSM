function [Population_Age,Pop_desire]=Age(Population_Age,Pop_desire,CasualProp)
% ageing all individuals who are alive by 1 day
x=0:1:50; %potential number of casual partners
ind_arr=find(Population_Age>0); % find all individuals who are alive
for ind=ind_arr
    Population_Age(ind)=Population_Age(ind)+1;

    if Population_Age(ind)<25
        ind0=1;
    elseif Population_Age(ind)<35
        ind0=2;
    elseif Population_Age(ind)<45
        ind0=3;
    elseif Population_Age(ind)<55
        ind0=4;
    elseif Population_Age(ind)<65
        ind0=5;
    else
        ind0=6;
    end
    cfX=cumsum(CasualProp(ind0,:));
    r2=rand;
    ind2=find(cfX>r2,1);
    Pop_desire(1,ind)=x(ind2);%desired number of partners within the last 12 months
end
end