function Pot_pairs=CreatePairList(candidates_ids)
% given a list of ids of individuals, this function creates a list of all
% possible pairs that can occur

num=size(candidates_ids,2);
Pot_pairs=zeros(num*(num-1)/2,2);

counter=1;
for counter1=1:1:(num-1)
    for counter2=(counter1+1):1:num
        Pot_pairs(counter,1)=candidates_ids(counter1);
        Pot_pairs(counter,2)=candidates_ids(counter2);
        counter=counter+1;
    end
end