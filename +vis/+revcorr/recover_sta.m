function sta = recover_sta(b)
%RECOVER_STA - recover the spike triggered average from the reconstruction block
%
% STA = vis.revcorr.recover_sta(B)
%
% Inputs:
%  B - the reconstruction block (MxMxTmax x N)
%
% Output:
%  STA - the spike triggered average (MxMxTmax)

block_sum = zeros(size(b,1), size(b,2), size(b,3));
for i = 1:size(b,4)
    block_sum = block_sum + b(:, :, :, i);
end
sta = block_sum / size(b,4);
end

