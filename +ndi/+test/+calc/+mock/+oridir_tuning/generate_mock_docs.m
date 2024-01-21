function output = generate_mock_docs(S)
% GENERATE_MOCK_DOCS - test generate mock documents for ndi.calc.vis.oridir_tuning
%
% OUTPUT = GENERATE_MOCK_DOCS(S)
%
% Generate mock docs for ndi.calc.vis.oridir_tuning calculator
% 

c = ndi.calc.vis.oridir_tuning(S);

scope = 'standard';

%[docs,doc_output,doc_expected_output] = c.generate_mock_docs('standard',9,'generate_expected_docs',1);
[docs,doc_output,doc_expected_output] = c.generate_mock_docs(scope,5);


b = [];
errormsg = {};
b_expected = [];

for i=1:numel(doc_output),
	for j=i:numel(doc_output),
		[doesitmatch,theerrormsg] = c.compare_mock_docs(doc_expected_output{i}, ...
			doc_output{i}, scope);
		b(i,j) = doesitmatch;
		b(j,i) = doesitmatch;
		errormsg{i,j} = theerrormsg;
		errormsg{j,i} = theerrormsg;
		[doesitmatch,theerrormsg] = c.compare_mock_docs(doc_expected_output{i},...
			doc_expected_output{j}, scope);
		b_expected(i,j) = doesitmatch;
		b_expected(j,i) = doesitmatch;
	end;
end;

output = vlt.data.var2struct('docs','doc_output','doc_expected_output','b','b_expected','errormsg');
