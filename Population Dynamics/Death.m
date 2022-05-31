function [Pop_Id,Pop_Age,Pop_Steady,Pop_history_steady,Rels_steady,Steady_dur,Pop_Casual,Pop_history_casual,Rels_casual,Casual_dur,Pop_desire]=Death(Pop_Id,Pop_Age,Pop_Steady,Pop_history_steady,DeathRates,Rels_steady,Steady_dur,Pop_Casual,Pop_history_casual,Rels_casual,Casual_dur,Pop_desire,tcounter)
    % retiring individuals from the population if the yearly mark has been
    % crossed and they have turned 83 or due to the background mortality

    % death due to the background mortality
    ind=find(Pop_Age>0);
    if numel(ind)>0
        age_backets=15:10:75;
        age_backets(end)=76;
        for counter=1:1:numel(ind)
            index=ind(counter);
            age=Pop_Age(index);
            age_bracket=find(age_backets>age,1)-1;
            mu=DeathRates(age_bracket);
            r1=rand;
            if r1<mu % individual dies
                id=Pop_Id(index);
                Pop_Id(index)=0;% remove the record from Id record
                Pop_Age(index)=0; % remove the record from age record
                % clean up steady relationships for the diseased and for their
                % partners
                SteadyParts=Pop_Steady(:,index);
                Pop_Steady(:,index)=0;
                Pop_history_steady(:,index)=-1;
                CasualPart=Pop_Casual(:,index);
                Pop_Casual(:,index)=0;
                Pop_history_casual(:,index)=-1;
                Pop_desire(:,index)=0;
                if SteadyParts(1)>0
                    ind_alter=find(Pop_Id==SteadyParts(1));
                    findego=find(Pop_Steady(:,ind_alter)==id);
                    Pop_Steady(findego,ind_alter)=0;
                    %clean up steady relationship list from this relationship
                    [tf,ind_rel]=ismember(sort([id,SteadyParts(1)]),Rels_steady(:,1:2),'rows');
                    Steady_dur=[Steady_dur,tcounter-Rels_steady(ind_rel,3)];
                    Rels_steady(ind_rel,:)=[];
                end
                % check the second possible partner
                if SteadyParts(2)>0
                    ind_alter=find(Pop_Id==SteadyParts(2));
                    findego=find(Pop_Steady(:,ind_alter)==id);
                    Pop_Steady(findego,ind_alter)=0;
                    %clean up steady relationship list from this relationship
                    [tf,ind_rel]=ismember(sort([id,SteadyParts(2)]),Rels_steady(:,1:2),'rows');
                    Steady_dur=[Steady_dur,tcounter-Rels_steady(ind_rel,3)];
                    Rels_steady(ind_rel,:)=[];
                end
                if sum(CasualPart>0)>0
                    ind_casual=find(CasualPart>0); % find all casual partners
                    for idalter=CasualPart(ind_casual) % go over all casual partners, one at a time
                        ind_alter=find(Pop_Id==idalter);
                        findego=find(Pop_Casual(:,ind_alter)==id);
                        Pop_Casual(findego,ind_alter)=0;
                        %clean up casual relationship list from this relationship
                        [tf,ind_rel]=ismember(sort([id,idalter]),Rels_casual(:,1:2),'rows');
                        Casual_dur=[Casual_dur,tcounter-Rels_casual(ind_rel,3)];
                        Rels_casual(ind_rel,:)=[];
                        Pop_desire(2,ind_alter)=Pop_desire(2,ind_alter)-1;
                    end
                end
            end
        end
    end

%     % death due to the being too old
    if mod(tcounter,365)==0 % a year has passed
        ind=find(Pop_Age>=75);
        if numel(ind)>0
            for counter=1:1:numel(ind)
                index=ind(counter);
                id=Pop_Id(index);
                Pop_Id(index)=0;% remove the record from Id record
                Pop_Age(index)=0; % remove the record from age record
                % clean up steady relationships for the diseased and for their
                % partners
                SteadyParts=Pop_Steady(:,index);
                Pop_Steady(:,index)=0;
                % clean up casual relationships for the diseased and their
                % partner
                CasualPart=Pop_Casual(:,index);
                Pop_Casual(:,index)=0;
                Pop_history_casual(:,index)=-1;
                Pop_desire(:,index)=0;
                if SteadyParts(1)>0
                    ind_alter=find(Pop_Id==SteadyParts(1));
                    findego=find(Pop_Steady(:,ind_alter)==id);
                    Pop_Steady(findego,ind_alter)=0;
                    %clean up steady relationship list from this relationship
                    [tf,ind_rel]=ismember(sort([id,SteadyParts(1)]),Rels_steady(:,1:2),'rows');
                    Steady_dur=[Steady_dur,tcounter-Rels_steady(ind_rel,3)];
                    Rels_steady(ind_rel,:)=[];
                end
                % check the second possible partner
                if SteadyParts(2)>0
                    ind_alter=find(Pop_Id==SteadyParts(2));
                    findego=find(Pop_Steady(:,ind_alter)==id);
                    Pop_Steady(findego,ind_alter)=0;
                    %clean up steady relationship list from this relationship
                    [tf,ind_rel]=ismember(sort([id,SteadyParts(2)]),Rels_steady(:,1:2),'rows');
                    Steady_dur=[Steady_dur,tcounter-Rels_steady(ind_rel,3)];
                    Rels_steady(ind_rel,:)=[];
                end
                if sum(CasualPart>0)>0
                    ind_casual=find(CasualPart>0); % find all casual partners
                    for idalter=(CasualPart(ind_casual))' % go over all casual partners, one at a time
                        ind_alter=find(Pop_Id==idalter);
                        findego=find(Pop_Casual(:,ind_alter)==id);
                        Pop_Casual(findego,ind_alter)=0;
                        %clean up steady relationship list from this relationship
                        [tf,ind_rel]=ismember(sort([id,idalter]),Rels_casual(:,1:2),'rows');
                        Casual_dur=[Casual_dur,tcounter-Rels_casual(ind_rel,3)];
                        Rels_casual(ind_rel,:)=[];
                        Pop_desire(2,ind_alter)=Pop_desire(2,ind_alter)-1;
                    end
                end
            end
        end
    end
end