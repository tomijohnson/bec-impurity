function iterate = checkconvergence(ait,acc)
    maxerror = max(abs(ait(end)-ait(end-3:end))/ait(end));
    if  maxerror < acc
        iterate = 0;
    else
        iterate = 1;
    end
end