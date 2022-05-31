function [suc_fl,Pop_casual,  Pop_casual_history, Pop_desire,Rels_casual]=create_casual_pairV2(t,Pop_Id,Pop_casual, Pop_casual_history,Rels_casual,Pop_desire,Pop_steady)
% The function updates tables People_casual for the two individuals, people history with the time stamp
% of initiation of the new partnership, register Rels with a new pair.
% ASSUMPTION ABOUT INPUT: ids in People_casual, People_history and
% Population_desire are in the same order
suc_fl=0;
r1=rand;
r2=rand;
indAlive=find(Pop_Id>0);
cs=cumsum(Pop_desire(1,indAlive))/sum(Pop_desire(1,indAlive));

ind1T=find(cs>r1,1);
ind2T=find(cs>r2,1);

ind1=indAlive(1,ind1T);
ind2=indAlive(1,ind2T);

id1=Pop_Id(1,ind1);
id2=Pop_Id(1,ind2);
% check that the same people or people already in a relationship were not
% picked
if id1~=id2 %not the same person
    parts1=Pop_casual(:,ind1);
    if ~ismember(id2,parts1) % not casual partners already
        parts1steady=Pop_steady(:,ind1);
        if ~ismember(id2,parts1steady)% id2 is not a steady partner of id1
            suc_fl=1;
            % find slots to insert new relationship
            insert_slot1=find(Pop_casual(:,ind1)==0,1);
            insert_slot2=find(Pop_casual(:,ind2)==0,1);
            Pop_casual(insert_slot1,ind1)=id2; % insert id of the new partner into the registry
            Pop_casual(insert_slot2,ind2)=id1; % insert id of the new partner into the new registry
    
            insert_slotT1=find(Pop_casual_history(:,ind1)==-1,1);
            insert_slotT2=find(Pop_casual_history(:,ind2)==-1,1);
            Pop_casual_history(insert_slotT1,ind1)=t;% record the time when the new relationship commenced for individual 1 
            Pop_casual_history(insert_slotT2,ind2)=t;% record the time when the new relationship commenced for individual 2
            
            Rels_casual=[Rels_casual;sort([id1, id2]),t]; % update relationship list Rels
            Pop_desire(2,ind1)=Pop_desire(2,ind1)+1; % increase the number of current partners
            Pop_desire(2,ind2)=Pop_desire(2,ind2)+1; % increase the number of the current partners
        end
    end
end
end