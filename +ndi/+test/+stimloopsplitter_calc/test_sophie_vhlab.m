function test_sophie_vhlab(S)

 % open a new experiment on the command line
 % S = ndi.session.dir([fullpathtoexperiment]); 

 % open a new object

c = ndi.calc.vis.stimloopsplitter_calc(S);

 % define our parameters
tic
parameters = c.default_search_for_input_parameters();

I = c.search_for_input_parameters(parameters);
toc %how much time does this take?
d = c.calculate(I{1});


d = c.calculate(I{2});


d = c.calculate(I{3});


d = c.calculate(I{4});



