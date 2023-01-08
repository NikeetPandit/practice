%Function finds aliased peak 
%Assuming only looking at positive frequency components
%Must know the sampled frequency and the true frequency (i.e., before is aliased)
%written by Nikeet Pandit

function f_aliased = find_aliased_peak(f_true, f_sample)

    f_nyquist = f_sample/2; 
    f_aliased_unfolded = f_true - f_sample;
    
    if f_true <= f_nyquist
        f_aliased = NaN;
    else
        f_aliased = mod(f_aliased_unfolded,f_nyquist); 
    end
    
    if f_true >= f_nyquist && f_aliased == 0
        f_aliased = -999; %Peak is missing in spectrum!
    end
    
end