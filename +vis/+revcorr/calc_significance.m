function [p_val, cdf, rescale] = calc_significance(sta,spike_num)
%CALC_SIGNIFICANCE - generate the p-value array and the rescaled matrix
%for plotting by taking in the sta and spike_num
%
% [P_VAL, RESCALE] = CALC_SIGNIFICANCE(STA,SPIKE_NUM)
%
% Inputs:
%  STA - the STA array for calculation
%  SPIKE_NUM, the number of spikes in the original json file
%
% OUTPUTS:
%  P_VAL - an array with the same size as STA array, and record the
%  significance of each pixel
%  RESCALE - an array with the same size as STA array. It takes the value
%  of P_VAL array and rescale it for plotting purposes. 

p_val = zeros(size(sta, 1), size(sta, 2), size(sta, 3));
cdf = zeros(size(sta, 1), size(sta, 2), size(sta, 3));
rescale = zeros(size(sta, 1), size(sta, 2), size(sta, 3));

for i = 1:size(sta,1)
	for j = 1:size(sta, 2)
		for k = 1:size(sta, 3)
			cdf(i,j,k) = vis.revcorr.get_cdf_val(sta(i,j,k), spike_num);
			if cdf(i,j,k) >= 0.5
				p_val(i,j,k) = -1 * log10(2 * (1 - cdf(i,j,k)));
			else
				p_val(i,j,k) = log10(2 * cdf(i,j,k));
			end
			if p_val(i,j,k) >= -log10(2*0.05),
				rescale(i, j, k) = vlt.math.rescale(p_val(i,j,k), [-log10(2*0.05) -log10(2*1e-8)], [193 256]);
			elseif p_val(i,j,k) <= log10(2*0.05)
				rescale(i, j, k) = vlt.math.rescale(p_val(i,j,k), [log10(2*1e-8) log10(2*0.05)], [0 64]);
			end
		end
	end
end

