function [left_pad,right_pad] = fix_pad(pad,gap_indx,a,w)
    % calculate left padding based on pad, a, w, and left-most index of gap
    if mod(gap_indx(1)-pad-1,w/a) == 0
        left_pad = pad;
    else
        left_pad = pad + mod(gap_indx(1)-pad-1,w/a);
    end
    % left_pad = pad;
    % left padding + gap_length + right padding must be a multiple of w/a
    rem = mod(left_pad+length(gap_indx)+left_pad,w/a);
    if rem == 0
        right_pad = left_pad;
    else
        if left_pad - rem < w/a
            right_pad = left_pad-rem+w/a;
        else
            right_pad = left_pad - rem;
        end
    end
end

