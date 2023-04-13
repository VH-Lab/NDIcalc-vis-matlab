function [b, errormsg] = oridir_compare_docs(document_expected, document_actual, scope)
% hartley_reverse_correlation_compare_docs
%
% [B, ERRORMSG] = ndi.calc.vis.test.hartley_reverse_correlation_compare_docs(DOC_EXPECTED, DOC_ACTUAL, SCOPE)
%


b = 1;
errormsg = '';

% start comparison

doc_e = document_expected.document_properties.hartley_reverse_correlation;
doc_a = document_actual.document_properties.hartley_reverse_correlation;

% Comparing stimulus_properties
% Comparing M

tol = %put in appropriate tolorance;

M_match = vlt.data.sizeeq(doc_e.stimulus_properties.M(:),doc_a.stimulus_properties.M(:));
if M_match
   M_match = max(abs(doc_e.stimulus_properties.M(:) - doc_a.stimulus_properties.M(:))) < tol;
end

if ~M_match
   b = 0;
   errormsg = ['M differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing L_max

tol = %put in appropriate tolorance;

L_max_match = vlt.data.sizeeq(doc_e.stimulus_properties.L_max(:),doc_a.stimulus_properties.L_max(:));
if L_max_match
   L_max_match = max(abs(doc_e.stimulus_properties.L_max(:) - doc_a.stimulus_properties.L_max(:))) < tol;
end

if ~L_max_match
   b = 0;
   errormsg = ['L max differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing K_max

tol = %put in appropriate tolorance;

K_max_match = vlt.data.sizeeq(doc_e.stimulus_properties.K_max(:),doc_a.stimulus_properties.K_max(:));
if K_max_match
   K_max_match = max(abs(doc_e.stimulus_properties.K_max(:) - doc_a.stimulus_properties.K_max(:))) < tol;
end

if ~K_max_match
   b = 0;
   errormsg = ['K max differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing sf_max

tol = %put in appropriate tolorance;

sf_max_match = vlt.data.sizeeq(doc_e.stimulus_properties.sf_max(:),doc_a.stimulus_properties.sf_max(:));
if sf_max_match
   sf_max_match = max(abs(doc_e.stimulus_properties.sf_max(:) - doc_a.stimulus_properties.sf_max(:))) < tol;
end

if ~sf_max_match
   b = 0;
   errormsg = ['sf max differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing fps

tol = %put in appropriate tolorance;

fps_match = vlt.data.sizeeq(doc_e.stimulus_properties.fps(:),doc_a.stimulus_properties.fps(:));
if fps_match
   fps_match = max(abs(doc_e.stimulus_properties.fps(:) - doc_a.stimulus_properties.fps(:))) < tol;
end

if ~fps_match
   b = 0;
   errormsg = ['fps differed by greater than ' num2str(tol) '.'];
   return;
end

%Comparing color_high

tol = %put in appropriate tolorance;

if any(abs(doc_e.stimulus_properties.color_high(:) - doc_a.stimulus_properties.color_high(:)) > tol)
   b = 0;
   errormsg = ['Color High differed by greater than ' num2str(tol) '.'];
end

%Comparing color_low

tol = %put in appropriate tolorance;

if any(abs(doc_e.stimulus_properties.color_low(:) - doc_a.stimulus_properties.color_low(:)) > tol)
   b = 0;
   errormsg = ['Color low differed by greater than ' num2str(tol) '.'];
end

%Comparing rect

tol = %put in appropriate tolorance;

if any(abs(doc_e.stimulus_properties.rect(:) - doc_a.stimulus_properties.rect(:)) > tol)
   b = 0;
   errormsg = ['Rect differed by greater than ' num2str(tol) '.'];
end

% Comparing reconstruction_properties
%Comparing T_coords

tol = %put in appropriate tolorance;

if any(abs(doc_e.reconstruction_properties.T_coords(:) - doc_a.reconstruction_properties.T_coords(:)) > tol)
   b = 0;
   errormsg = ['T_coords differed by greater than ' num2str(tol) '.'];
end

%Comparing X_coords

tol = %put in appropriate tolorance;

if any(abs(doc_e.reconstruction_properties.X_coords(:) - doc_a.reconstruction_properties.X_coords(:)) > tol)
   b = 0;
   errormsg = ['X_coords differed by greater than ' num2str(tol) '.'];
end

%Comparing Y_coords

tol = %put in appropriate tolorance;

if any(abs(doc_e.reconstruction_properties.Y_coords(:) - doc_a.reconstruction_properties.Y_coords(:)) > tol)
   b = 0;
   errormsg = ['Y_coords differed by greater than ' num2str(tol) '.'];
end

% Going into superclass
% Comparing reverse_correlation
% Comparing method

method_match = strcmpi(doc_e.superclass.reverse_correlation.method, doc_a.superclass.reverse_correlation.method);
if ~method_match
   b = 0;
   errormsg = ['Expected superclass.reverse_correlation.method field of ' doc_e.superclass.reverse_correlation.method ' but observed ' doc_a.superclass.reverse_correlation.method];
   return;
end

% Comparing dimension_labels

if any(doc_e.superclass.reverse_correlation.dimension_labels(:) ~= doc_a.fit.superclass.reverse_correlation.dimension_labels(:))
   b = 0;
   errormsg = ['Expected superclass.reverse_correlation.dimension_labels field of ' doc_e.superclass.reverse_correlation.dimension_labels ' but observed ' doc_a.superclass.reverse_correlation.dimension_labels];
   return;
end

end
