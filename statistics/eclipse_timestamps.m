timestamps= [];
e =1;
for i = 1:length(eclipse)
    if eclipse(i) == e
        timestamps = [timestamps i];
        if e == 1
            e = 0;
        else
            e = 1;
        end
    end
end
timestamps