function test_sophie_vhlab(S)

 % open a new experiment on the command line
 % S = ndi.session.dir([fullpathtoexperiment]); 

 % open a new object
c = ndi.calc.vis.stimloopsplitter(S);

 % define our parameters

parameters = c.default_search_for_input_parameters();

I = c.search_for_input_parameters(parameters);

d = c.calculate(I{1});


d = c.calculate(I{2});


d = c.calculate(I{3});


d = c.calculate(I{4});



