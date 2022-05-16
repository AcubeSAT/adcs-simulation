function accuracy_error(R_k,cov_threshold)
    for i = 1:length(R_k)
        if R_k(i) > cov_threshold(i)
            faulty_measurement = i;
            ME = MException('ErrorFound:AccuracyErrorOccured','error on sensor %d',faulty_measurement);
            throw(ME);
        end
    end
    
end