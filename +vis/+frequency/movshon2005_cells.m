function [cell1, cell2] = movshon2005_cells()
% MOVSHON2005_CELLS - reproduced data from 2 examples cells from Movshon et al
%
% [CELL1,CELL2] = MOVSHON2005_CELLS()
%
% Return temporal frequency vs. response rates for 2 example cells
% from Figure 1a from Movhson et al. 2005 (J Neurosci., V 25).
% 
% Data grabbed with DataThief.
%
% Example:
%   [cell1,cell2] = vis.frequency.movshon2005_cells();
%   figure;
%   plot(cell1(:,1),cell1(:,2),'ro-');
%   hold on;
%   plot(cell2(:,1),cell2(:,2),'bo-');
%   set(gca,'xscale','log');
%   xlabel('Temporal frequency (Hz)');
%   ylabel('Response (ips)');
%   box off;
%   
% 


cell1 = [ 0.3243, 11.7454
0.6524, 16.0866
1.3415, 19.9115
2.5354, 15.974
5.0205, 0.6593];

[dummy,cell1_indexes] = sort(cell1(:,1));

cell1 = cell1(cell1_indexes,:);

cell2 = [ ...
1.7191, 12.3632
3.4365, 14.2786
6.7911, 4.5274
13.5011, 0.8706
27.234, 0.2612
0.428, 3.3955 ];

[dummy,cell2_indexes] = sort(cell2(:,1));

cell2 = cell2(cell2_indexes,:);
