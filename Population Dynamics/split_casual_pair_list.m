function [Pop_casual, Rels_casual,Rels_dur,Pop_desire]=split_casual_pair_list(Pop_Id,Pop_casual,Rels_casual,Rels_dur,Pairs_dissolve,Pop_desire,tcounter)
%dissolve the pairs whose ids are provided in Pairs_dissolve list, update
%tables People_casual,People_history,Rels,Population_desire,Rels_dur
 Rels_dissolve=Rels_casual(Pairs_dissolve,:);
 Rels_casual(Pairs_dissolve,:)=[]; % remove the pairs from the registry
 % calculate the durations of the relationships that are about to be
 % dissolved
 durations=tcounter-Rels_dissolve(:,3);
 Rels_dur=[Rels_dur;durations]; %record the durations of terminated relationships
 for counter=1:1:size(Rels_dissolve,1)
    % retrieve ids of the couple
    id1=Rels_dissolve(counter,1);
    id2=Rels_dissolve(counter,2);
    % remove alters from the partner list for the two egos
    ind1=find(Pop_Id==id1); % find the location of id1 data in all arrays
    ind2=find(Pop_Id==id2); % find the location of id2 data in all arrays
    id1alter=find(Pop_casual(:,ind1)==id2);
    id2alter=find(Pop_casual(:,ind2)==id1);
    Pop_casual(id1alter,ind1)=0;
    Pop_casual(id2alter,ind2)=0;
    Pop_desire(2,ind1)=Pop_desire(2,ind1)-1;
    Pop_desire(2,ind2)=Pop_desire(2,ind2)-1;
 end
end