function [low,maxv,high]= compute_halfwidth( x, y )

% COMPUTE_HALFWIDTH - find where signal reaches maximum, 1/2 on left side, 1/2 on right side
%
%     [LOW,MAX,HIGH] = COMPUTE_HALFWIDTH( X, Y )
%
%     returns MAX, where x position where function attains its maximum value
%     LOW < MAX,  where function attains half its maximum
%     HIGH > MAX, where function attains half its maximum
%     returns NAN for LOW or/and HIGH, when function does not come below the point
%
%     See also: COMPUTE_HALFWIDTH_INTERP

[maxvalue,max_loc] = max(y);
halfheight = maxvalue/2;
maxv = x(max_loc);

before = y(max_loc-1:-1:1);
after = y(max_loc+1:end);

low_loc = find(before<=halfheight,1);
high_loc = find(after<=halfheight,1);

if isempty(low_loc),
	low = NaN;
else,
	if (max_loc-low_loc)<= 0, keyboard; end;
	low = x( max_loc-low_loc );
end;

if isempty(high_loc),
	high = NaN;
else,
	high = x( max_loc+high_loc );
end;

