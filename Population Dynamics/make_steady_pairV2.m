function [s_fl,Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady]=make_steady_pairV2(Pop_Id,Pop_Steady,Pop_steady_history,Rels_steady,psi,Pop_casual,tcounter)
% this function creates a pair such that individuals with no existing
% steady pairs have 1/psi probability times chance higher of being picked
% input: Population_Id [id1,id2...] Population_Steady table [part1_id;part2_id] detailing steady partners
% of each individual, no more than 2 are possible. Casual_Population
% population table [id,part1_id,part2_id,...,partMax_casual_id] detailing
% casual parters of each individual, no more than Max_casual casual
% partners is possible. Rels - relationship table capturing all
% relationships in the population

%Eventually will need to return Population_casual and Max_casual

% check of the input
ownName=mfilename;

if psi<0 | psi>1
    error([ownName,': psi should be between 0 and 1']);
end

s_fl=0;

%establish existing individuals and their locations
ind=find(Pop_Id>0);
% calculate the number of partners for each person

num_parts=(Pop_Steady(1,ind)>0)+(Pop_Steady(2,ind)>0);
% set propensities to be chosen
a=(num_parts==0)+psi*(num_parts==1)+0*(num_parts==2);
% normalize
if sum(a)>0 % there is someone to pick
    a=a/sum(a);
    cs=cumsum(a);
    r1=rand;
    r2=rand;

    indT1=find(cs>r1,1);
    indT2=find(cs>r2,1);
    
    if indT1~=indT2 % the person was not selected twice
        %find indices of the chosen individuals
        ind_id1=ind(indT1);
        ind_id2=ind(indT2);
        % find the id of the first person
        id1=Pop_Id(1,ind_id1);
        % find the id of the second person
        id2=Pop_Id(1,ind_id2);
        partners_id1=Pop_Steady(1:2,ind_id1);% partners of id1
        partners_id2=Pop_Steady(1:2,ind_id2);% partners of id2
        num_pars1=2-numel(find(partners_id1==0));
        num_pars2=2-numel(find(partners_id2==0));

        % check for casual partners
        partners_casual_id1=Pop_casual(:,ind_id1);
        if sum(partners_casual_id1==id2)==0
            if num_pars1==0 & num_pars2==0 % both people do not have partners
                % create the pair
                Pop_Steady(1,ind_id1)=id2;
                Pop_Steady(1,ind_id2)=id1;
                Pop_steady_history(1,ind_id1)=tcounter;
                Pop_steady_history(1,ind_id2)=tcounter;
                Rels_steady=[Rels_steady; sort([id1,id2]) tcounter];
                s_fl=1;
            elseif num_pars1<2 & num_pars2<2 % both individuals can have another steady relationship but at least one of them already has a partner
                if num_pars1==1 & num_pars2==0
                    % find for the individual who already has a
                    % relationship an available slot for the new partner
                    ind_t=find(Pop_Steady(:,ind_id1)==0);
                    Pop_Steady(ind_t,ind_id1)=id2;
                    Pop_Steady(1,ind_id2)=id1;
                    Pop_steady_history(ind_t,ind_id1)=tcounter;
                    Pop_steady_history(1,ind_id2)=tcounter;
                    Rels_steady=[Rels_steady; sort([id1,id2]) tcounter];
                    s_fl=1;
                elseif num_pars1==0 & num_pars2==1
                    % find for the individual who already has a
                    % relationship an available slot for the new partner
                    ind_t=find(Pop_Steady(:,ind_id2)==0);
                    Pop_Steady(1,ind_id1)=id2;
                    Pop_Steady(ind_t,ind_id2)=id1;
                    Pop_steady_history(1,ind_id1)=tcounter;
                    Pop_steady_history(ind_t,ind_id2)=tcounter;
                    Rels_steady=[Rels_steady; sort([id1,id2]) tcounter];
                    s_fl=1;
                else
                    if sum(Pop_Steady(1:2,ind_id1)==id2)==0
                        % find for both individuals an available slot for the new partner
                        ind_t1=find(Pop_Steady(:,ind_id1)==-1);
                        ind_t2=find(Pop_Steady(:,ind_id2)==-1);
                        Pop_Steady(ind_t1,ind_id1)=id2;
                        Pop_Steady(ind_t2,ind_id2)=id1;
                        Pop_steady_history(ind_t1,ind_id1)=tcounter;
                        Pop_steady_history(ind_t2,ind_id2)=tcounter;
                        Rels_steady=[Rels_steady; sort([id1,id2]) tcounter];
                        s_fl=1;
                    end
                end
            end
        end

    end
end