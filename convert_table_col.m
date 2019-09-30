% Script to convert columns of tables in workspace 
%   from char to double (Xiaoyi's CF files)

column='year';


% get all the variables
list=who;

% loop over variable list
for i=1:length(list)
    
    % check if variable is a table
    if istable(eval(list{i}))
        
        % check if table has specified column
        if any(strcmp(column,eval([list{i} '.Properties.VariableNames'])))
            
            % convert column from char to double
            eval([list{i} '.' column '=str2num(' list{i} '.' column ');']);
        end
    end
end

clearvars column list i 