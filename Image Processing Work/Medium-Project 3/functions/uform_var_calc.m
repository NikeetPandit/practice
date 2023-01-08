%--- Function to calculate uniform noise variance
%
%Written by Nikeet Pandit

function var_uform = uform_var_calc(a,b)
    var_uform = ((b-a)^2)/12; 
end