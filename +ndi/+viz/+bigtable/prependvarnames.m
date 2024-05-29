function T_out = prependvarnames(T,prefix)
% PREPENDVARNAMES - add a prefix to variable names of a table
% 
% T_OUT = PREPENDVARNAMES(T, PREFIX)
% 

v = prefix + T.Properties.VariableNames;
T_out = renamevars(T,T.Properties.VariableNames,v);
