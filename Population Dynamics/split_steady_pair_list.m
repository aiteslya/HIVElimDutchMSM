function [Population_Steady,Rels_Steady,Steady_dur]=split_steady_pair_list(Population_Id,Population_Steady,Rels_Steady,pairs_list,tcounter,Steady_dur)
% function splits pairs whose ids are in the pair list, removes the pair from the
% relationship table Rels and updates the table People to remove the
% relationship for the respective individuals

%retrieve pairs that need to be split
RelsToSplit=Rels_Steady(pairs_list,:);
num_pairs=size(RelsToSplit,1);
% remove from rels table
Rels_Steady(pairs_list,:)=[];


for counter=1:1:num_pairs
    Steady_dur=[Steady_dur, tcounter-RelsToSplit(counter,3)];
    %individuals who are splitting
    id1=RelsToSplit(counter,1);
    id2=RelsToSplit(counter,2);
    
    % remove id2 from the registry for id1
    indego1=find(Population_Id(1,:)==id1);
    indalter1=find(Population_Steady(1:2,indego1)==id2);
    if indalter1==1
        if Population_Steady(2,indego1)==0 % no other partners
            Population_Steady(indalter1,indego1)=0;
        else % one more partner, we need to shift it up
            tempalter=Population_Steady(2,indego1);
            Population_Steady(1,indego1)=tempalter;
            Population_Steady(2,indego1)=0;
        end    
    else    
        Population_Steady(indalter1,indego1)=0;
    end
    
    % remove id1 from the registry for id2
    indego2=find(Population_Id(1,:)==id2);
    indalter2=find(Population_Steady(1:2,indego2)==id1);
    if indalter2==1
        if Population_Steady(2,indego2)==0 % no other partners
            Population_Steady(indalter2,indego2)=0;
        else % one more partner, we need to shift it up
            tempalter=Population_Steady(2,indego2);
            Population_Steady(1,indego2)=tempalter;
            Population_Steady(2,indego2)=0;
        end    
    else    
        Population_Steady(indalter2,indego2)=0;
    end

end
end

