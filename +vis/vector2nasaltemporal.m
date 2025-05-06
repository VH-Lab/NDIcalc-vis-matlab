function v_out = vector2nasaltemporal(v_in, viewing_eye)
% VECTOR2NASALTEMPORAL - convert a vector from eye view to up/down/nasal/temporal view
%
% V_OUT = ndi.fun.vis.vector2nasaltemporal(V_IN, VIEWING_EYE)
%
% Given a vector V_IN in degrees in compass coordinates with VIEWING_EYE
% either 'left' or 'right, converts the vector to a coordinate system
% such that the left half of the viewing field corresponds to 'temporal' and 
% the right half of the viewing field corresponds to 'nasal'.
%
% Inputs: V_IN (degrees in compass coordinates)
%         VIEWING_EYE (can be 0 or false (left) or 1 or true (right), or the strings 'left' or 'right')
% 
% Example:
%   v_in = 90; % a right-ward vector
%   viewing_eye = 'right'; % this vector moves temporally
%   v_out = ndi.fun.vis.vector2nasaltemporal(v_in,viewing_eye);
%      % v_out == 270
% 

deg_cartesian = vlt.math.compass2cartesian(v_in,0);
rad_cartesian = vlt.math.deg2rad(deg_cartesian);
[x,y] = pol2cart(rad_cartesian,1);

if ischar(viewing_eye),
	if isempty(find(strcmpi(viewing_eye,{'left','right'}))),
		error(['VIEWING_EYE must be LEFT or RIGHT when passed as a character string.']);
	end;
else,
	viewing_eye = double(viewing_eye); % in case of boolean inputs
end;

if isnumeric(viewing_eye),
	if viewing_eye==0,
		viewing_eye = 'left';
	elseif viewing_eye==1,
		viewing_eye = 'right';
	else,
		error(['Unknown numeric input for VIEWING_EYE.']);
	end;
end;

v_out = [];

if strcmpi(viewing_eye,'left'),
	v_out = v_in;
elseif strcmpi(viewing_eye,'right'),
	x = -x;
	[th,r] = cart2pol(x,y);
	v_out_cart = vlt.math.rad2deg(th);
	v_out = vlt.math.cartesian2compass(v_out_cart,0);
else,	
	error(['Unknown VIEWING_EYE input ' viewing_eye '.']);
end;

v_out = mod(v_out,360);
v_out = mod(v_out,360); % there's a possibility the above can produce 360 as an output , if input is -1e-14 for example)

