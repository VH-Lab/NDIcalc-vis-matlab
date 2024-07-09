function p = speed_nested_f(response_vector_length,error_withspeed,error_nospeed)
% SPEED.FIT.SPEED_NESTED_F  Perform F-test on models varying in parameters.
%
%  P = SPEED.FIT.SPEED_NESTED_F(RESPONSE_VECTOR_LENGTH,ERROR_WITHSPEED,ERROR_NOSPEED)
%
%  Performs an F-test between two models with and without a speed parameter
%  to determine whether a neuron exhibits significant speed tuning.
%
%  Inputs:
%    RESPONSE_VECTOR_LENGTH is the value of the lenght of the vector of measured responses driven by the spatial and
%    temporal frequency values.
%
%    ERROR_WITHSPEED is the value of the squared 2-norm of the residual
%    of the model with a speed parameter.
%
%    ERROR_NOSPEED is the value of the squared 2-norm of the residual
%    of the model without a speed parameter.
%
%  Outputs:
%    P is the p-value returned by the F-test performed between the two
%    models.
%
%  See also SPEED.FIT.FIT, SPEED.FIT.FIT_NOSPEED.
%
    df_withspeed = response_vector_length - 7; % The number of degrees of freedom in the model with a speed parameter
    df_nospeed = response_vector_length - 6; % The number of degrees of freedom in the model without a speed parameter
    SSE_withspeed = error_withspeed; % The error sum of squares of the model with a speed parameter
    SSE_nospeed = error_nospeed; % The error sum of squares of the model without a speed parameter

    changeinerror = (SSE_nospeed - SSE_withspeed) / SSE_withspeed;
    changeindegreesoffreedom = (df_nospeed - df_withspeed) / df_withspeed;

    F = changeinerror / changeindegreesoffreedom; % The F-statistic
    p = 1 - fcdf(F,df_nospeed-df_withspeed,df_withspeed); % The p-value from the F-test
