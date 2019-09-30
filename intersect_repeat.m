function [ C, i_x1, i_x2 ] = intersect_repeat( x1,x2 )
%INTERSECT_REPEAT like intersect, but keeps repetitions
%
% Input:
%   x1,x2: two arrays with the same number of columns
%
% Output: 
%   C: rows common to x1 and x2
%           ---
%           If a row repeats in one array, C contains that row the
%           same number of times it was repeated
%           ---
%           If a row repeats in both arrays, C contains that row a
%           number of times equal to the max n.o. permutations (repeats in
%           x1 times repeats in x2) -- intended for matching larger arrays,
%           where repetitions in x1,x2 still correspond to unique lines in
%           the original array
%
%   i_x1, i_x2: indices such that x1(i_x1,:) = x2(i_x2,:) = C
%
%
% Kristof Bognar, February 2018



%% find matches

% matches in array 1
[~,ind12]=ismember(x1,x2,'rows');

% matches in array 2
[~,ind21]=ismember(x2,x1,'rows');

% stop if there are no matches
if sum(ind12)==0
   
    C=[]; i_x1=[]; i_x2=[];
    
    return
    
end


%% unique matches and first occurence of each element

[tmp1,ia1,~]=unique(ind12);
[tmp2,~,~]=unique(ind21);

% find how many times each index repeats
reps_unique1=hist(ind12,tmp1);
reps_unique2=hist(ind21,tmp2);

% remove zeros (no match)
ia1(tmp1==0)=[];
reps_unique1(tmp1==0)=[];
tmp1(tmp1==0)=[];

reps_unique2(tmp2==0)=[];
tmp2(tmp2==0)=[];


%% loop over matches
inds_final=[0,0];

% check elements from one side only, the problem is symmetric
for i=1:length(ia1)

    % number of matches in x2 for current x1 line
    x2reps=reps_unique2(tmp2==ia1(i));
    
    % if matching element only appears once in x1
    if reps_unique1(i)==1
        
        if x2reps==1
            % one to one match
            inds_final=[inds_final; [ia1(i), tmp1(i)]];
            
        else
            % x2 contains the repetition
            % find indices of repeated value in x2
            ind_reps=find(ind21==ia1(i));
            inds_final=[inds_final; [ia1(i)*ones(size(ind_reps)), ind_reps] ];
            
        end
    
    else % if matching element repeats in x1
        
        if x2reps==1
            % x1 contains the repetition
            % find indices of repeated value in x1
            ind_reps=find(ind12==tmp1(i));
            inds_final=[inds_final; [ind_reps, tmp1(i)*ones(size(ind_reps))] ];
            
        else
            % both arrays have repeated values
            ind_reps2=find(ind21==ia1(i));
            ind_reps1=find(ind12==tmp1(i));
            
            tmp=length(ind_reps1)*length(ind_reps2);
            
            inds_final=[ inds_final;...
                         [repelem(ind_reps1, tmp/length(ind_reps1)),...
                          repmat(ind_reps2, tmp/length(ind_reps2), 1)]    ];
            
        end
        
    end

end

% remove dummy first line
inds_final(1,:)=[];

% sort indices
inds_final=sortrows(inds_final);

%% final output

i_x1=inds_final(:,1);
i_x2=inds_final(:,2);

C=x1(i_x1,:);

% % final check
% if ~isequal(C,x2(i_x2,:))
%     error('Something went terribly wrong')
% end

end

