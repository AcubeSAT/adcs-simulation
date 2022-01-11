real_v = importdata('real_v.csv');
bdot_v = importdata('bdot_v.csv');
D2N_threshold = 0.0035;

counter = 0;
i = 0;
while counter < 50
    i = i + 1;
    if (sum(real_v(:,i) < D2N_threshold) == 3)
        counter = counter + 1;
    else
        counter = 0;
    end
end
true_flag = i;

N = 400;



%% CASE 1 - Mean of Velocity
flag_bdot_mean = 0;
index = 0;

while (flag_bdot_mean == 0)
    index = index + 1;
    
    matrix_with_real_v = real_v(:,index:index + N - 1);
    matrix_with_bdot_v = bdot_v(:,index:index + N - 1);
    
    mean_real_v = mean(matrix_with_real_v, 2);
    mean_bdot_v = mean(matrix_with_bdot_v, 2);

    if (sum(mean_real_v < D2N_threshold) == 3)
        flag_real_mean = 1;
    else
        flag_real_mean = 0;
    end

    if (sum(mean_bdot_v < D2N_threshold) == 3)
        flag_bdot_mean = 1;
    else
        flag_bdot_mean = 0;
    end

end

flag_index_mean = index;

%% CASE 2 - Median of Velocity

flag_bdot_median = 0;
index = 0;

while (flag_bdot_median == 0)
    
    index = index + 1;
    
    matrix_with_real_v = real_v(:,index:index + N - 1);
    matrix_with_bdot_v = bdot_v(:,index:index + N - 1);
    
    median_real_v = median(matrix_with_real_v, 2);
    median_bdot_v = median(matrix_with_bdot_v, 2);

    if (sum(median_real_v < D2N_threshold) == 3)
        flag_real_median = 1;
    else
        flag_real_median = 0;
    end

    if (sum(median_bdot_v < D2N_threshold) == 3)
        flag_bdot_median = 1;
    else
        flag_bdot_median = 0;
    end
end

flag_index_median = index;

%% CASE 3 - Trigger processing the length of the vector

limit = 350;
flag_bdot_no_exceptions = 0;

while (flag_bdot_no_exceptions == 0)
    
    index = index + 1;
    
    matrix_with_real_v = real_v(:,index:index + N - 1);
    matrix_with_bdot_v = bdot_v(:,index:index + N - 1);

    compare_real_v = matrix_with_real_v < D2N_threshold;
    compare_bdot_v = matrix_with_bdot_v < D2N_threshold;

    for i = 1:3
        if (sum(compare_real_v(i,:)) < limit)
            flag_real_no_exceptions = 0;
            break
        else
            flag_real_no_exceptions = 1;
        end
        if (sum(compare_bdot_v(i,:)) < limit)
            flag_bdot_no_exceptions = 0;
            break
        else
            flag_bdot_no_exceptions = 1;
        end
    end
end

flag_index_length = index;

%% CASE 4 - Trigger with exceptions on how many consecutive steps can be different

limit = 350;
exceptions = 10;
flag_bdot_exceptions = 0;

while (flag_bdot_exceptions == 0) && index < N
    
    index = index + 1;
    
    matrix_with_real_v = real_v(:,index:index + N - 1);
    matrix_with_bdot_v = bdot_v(:,index:index + N - 1);

    compare_real_v = matrix_with_real_v < D2N_threshold;
    compare_bdot_v = matrix_with_bdot_v < D2N_threshold;

    for i = 1:3
        if (sum(compare_real_v(i,:)) < limit)
            flag_real_exceptions = 0;
            break
        else
            for j = 1:exceptions:N-exceptions
                if (sum(compare_real_v(i,j:j+exceptions) == 0))
                    flag_real_exceptions = 0;
                    break
                else
                    flag_real_exceptions = 1;
                end
            end
        end
    end

    for i = 1:3
        if (sum(compare_bdot_v(i,:)) < limit)
            flag_bdot_exceptions = 0;
            break
        else
            for j = 1:exceptions:N-exceptions
                if (sum(compare_bdot_v(i,j:j+exceptions) == 0))
                    flag_bdot_exceptions = 0;
                    break
                else
                    flag_bdot_exceptions = 1;
                end
            end
        end
    end
end

flag_index_exceptions = index;