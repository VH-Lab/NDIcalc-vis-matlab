function [b,msg] = values_within_tolerance(v1, v2, tolerance, fieldname)
% VALUES_WITHIN_TOLERANCE - check to see if values are within a tolerance
% 
% [B,MSG] = VALUES_WITHIN_TOLERANCE(V1, V2, TOLERANCE, FIELDNAME)
%
% This function performs a number of checks on the values of arrays V1 and V2.
%  a) It checks to make sure the arrays are the same size.
%  b) It checks to make sure no values differ by an amount greater than TOLERANCE 
%
% If both tests pass, then B is 1 and MSG is blank/empty ('').
% If either test fails, then B is 0 and MSG contains an error message as to 
% how it failed.
% 
% Example:
%  v1 = [ 1 2 3 ];
%  v2 = [ 1 2 3 ];
%  v3 = [ 1 2 ];
%  v4 = [ 1 2 4];
%
% [b1,msg1] = ndi.test.values_within_tolerance(v1, v2, 0, 'myfield')
%  % b1 is 1, msg1 is ''
% [b2,msg2] = ndi.test.values_within_tolerance(v1, v3, 0, 'myfield')
%  % b2 is 0, msg2 is 'Arrays of myfield are not the same size'
% [b3,msg3] = ndi.test.values_within_tolerance(v1, v4, 0, 'myfield')
%  % b3 is 0, msg3 is 'Differences in arrays of myfield exceed the tolerance provided (0)'];

b = 1;
msg = '';

size_matched = vlt.data.sizeeq(v1,v2);

if ~size_matched,
  b = 0;
  msg = ['Arrays of ' fieldname ' are not the same size.'];
  return;
end;

if isnan(v1)
    b = 0;
    msg = [fieldname ' is not a number'];
    return
end

if isnan(v2)
    b = 0;
    msg = [fieldname ' is not a number'];
    return
end

tolerance_matched = max(abs(v1(:) - v2(:))) < tolerance;

if ~tolerance_matched,
   b = 0;
   msg = ['Differences in arrays of ' fieldname ' exceed the tolerance provided ' ( num2str(tolerance))];
end;
   
end
