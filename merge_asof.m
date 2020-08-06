function [ table_out ] = merge_asof( table1, table2, merge_on, col_t2, dt, col_append )
%MERGE_ASOF( table1, table2, merge_on, col_t2, dt )
%
% after the pandas function in python -- left join on nearest value
% for each value in table1, finds the nearest line in table 2, and adds the
% value corresponding value from the selected column
%
% INPUT:
%   table1, table2: tables with any columns
%   merge_on: column to use when finding nearest values (e.g. fractional
%             time); must be present in both tables
%   col_t2: column in table2 that is to be added to table1
%   dt: tolerance for the nearest values (units must match the merge_on column)
%       if there are no matches within the tolerance, the meged value is NaN
%   col_append: the new column in table1 will be called [col_t2 col_append]
%
%@Kristof Bognar, August 2020

% find differences between each element
% bsxfun uses less memory (and requires fewer lines of code) compared to repmat
diff=abs(bsxfun(@minus,table2.(merge_on),table1.(merge_on)'));
% each column represents the difference of one table1 time to all the table2 times
% each row represents the difference of one table2 time to all the table1 times

% indices of minima in table2, corresponding to each table1 element
% (closest times to each table1 time)
[min_val12,min_ind12]=min(diff,[],1);

% output is table 1 with extra column
table_out=table1;
table_out.([col_t2 col_append])=table2.(col_t2)(min_ind12);

% remove merged values outside of time window
table_out.([col_t2 col_append])(min_val12>dt)=NaN;

end

