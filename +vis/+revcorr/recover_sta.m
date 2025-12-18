function sta = recover_sta(b)
%RECOVER_STA Summary of this function goes here
%   Detailed explanation goes here
block_sum = zeros(size(b,1), size(b,2), size(b,3));
for i = 1:size(b,4)
    block_sum = block_sum + b(:, :, :, i);
end
sta = block_sum / size(b,4);
end

