% function [out1, out2, ] = function_name(in1, in2)

function [max_val, max_loc] = my_max(v)
    n = length(v);
    max_val = -inf;
    
    for i = 1:n
        if v(i) > max_val
            max_val = v(i);
            max_loc = i;
        end
    end