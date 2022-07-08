function docs = oridir_generate_mock_docs(type, number)
% ndi.calc.vis.test.oridir_generate_mock_docs - generate mock input documents for testing ndi.calc.vis.oridir_tuning
%
% DOCS = ndi.calc.vis.test.oridir_generate_mock_docs(TYPE, NUMBER)
%
% Generate mock input documents depending upon the type called for. TYPE can be 'standard', 'random_nonoise', or 'random_noisy'.
% NUMBER specifies the number of documents to generate
switch type
    case standard
        for x = 1.0:+1.0:number
            P(x,:)= [x mod(x*10.0,40) 0.5*mod(x*10.0,40) 90-mod(x,2)*45 x*30]; %This is the standard preset value
        end
   
    case random_noisy
        for x = 1.0:+1.0:number
            P(x,:)= ndi.calc.vis.test.oridir_mock_parameters; %This is the parameters with noise, I used the model 10% of response + 2
            P(x,1) = P(x,1)*1.1 + 2;
            P(x,2) = P(x,2)*1.1 + 2;
            P(x,3) = P(x,3)*1.1 + 2;
        end

    case random_nonoise
        for x = 1.0:+1.0:number
            P(x,:)= ndi.calc.vis.test.oridir_mock_parameters; %This is parameters without noise pulled straight from the function ndi.calc.vis.test.oridir_mock_parameters
        end

end
