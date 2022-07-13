function docs = oridir_generate_mock_docs(type, number)
% ndi.calc.vis.test.oridir_generate_mock_docs - generate mock input documents for testing ndi.calc.vis.oridir_tuning
%
% DOCS = ndi.calc.vis.test.oridir_generate_mock_docs(TYPE, NUMBER)
%
% Generate mock input documents depending upon the type called for. TYPE can be 'standard', 'random_nonoise', or 'random_noisy'.
% NUMBER specifies the number of documents to generate
switch type
    case standard
        % standard should be a set list, it can be sampled over and over again to make a list more than 10
        % parameters are [Rsp Rp Rn sig op] as described in Mazurek et al.
        P_(1,:) = [ 0 20 10 30 45] ; % response of 20 in preferred direction of 45 degrees, 10 opposite
        P_(2,:) = [ 0 20 10 45 45] ; % broader tuning
        P_(3,:) = [ 0 20 10 90 45] ; % really broad tuning 
        P_(4,:) = [ 0 20 10 90 45] ; % really broad tuning 
        P_(5,:) = [ 10 20 10 30 45] ; % large offset
        P_(6,:) = [ 10 20 19 30 45] ; % really low direction index offset
        % try to increase to 10-15 cases, using Mazurek et al. figure to pick good variation
        for x = 1.0:+1.0:number
            index = 1 + mod(x-1,size(P_,2));
            P(x,:)= P_(index,:);
        end
   
    case random_noisy
        for x = 1.0:+1.0:number
            P(x,:)= ndi.calc.vis.test.oridir_mock_parameters; % This is the parameters with noise, I used the model 10% of response + 2. I'm not quite sure is this the 
            P(x,1) = P(x,1)*1.1 + 2;                            % right way to call a function in MATLAB
            P(x,2) = P(x,2)*1.1 + 2;
            P(x,3) = P(x,3)*1.1 + 2;
        end

    case random_nonoise
        for x = 1.0:+1.0:number
            P(x,:)= ndi.calc.vis.test.oridir_mock_parameters; % This is parameters without noise pulled straight from the function ndi.calc.vis.test.oridir_mock_parameters
        end

end
