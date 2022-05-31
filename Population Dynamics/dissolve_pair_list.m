function [People_casual, Rels,Rels_dur]=dissolve_pair_list(t,People_casual,Rels,Rels_dur,Pairs_dissolve)
%dissolve the pairs whose ids are provided in Pairs_dissolve list, update
%tables People_casual,People_history,Rels,Population_desire,Rels_dur
 Rels_dissolve=Rels(Pairs_dissolve,:);
 Rels(Pairs_dissolve,:)=[]; % remove the pairs from the registry
 % calculate the durations of the relationships that are about to be
 % dissolved
 durations=t-Rels_dissolve(:,3);
 Rels_dur=[Rels_dur;durations]; %record the durations of terminated relationships
 for counter=1:1:size(Rels_dissolve,1)
    % retrieve ids of the couple
    id1=Rels_dissolve(counter,1);
    id2=Rels_dissolve(counter,2);
    % remove alters from the partner list for the two egos
    ind1=find(People_casual(1,:)==id1);
    ind2=find(People_casual(1,:)==id2);
    id1alter=find(People_casual(2:end,ind1)==ind2)+1;
    id2alter=find(People_casual(2:end,ind2)==ind1)+1;
    People_casual(id1alter,ind1)=0;
    People_casual(id2alter,ind2)=0;
 end
end