
    function forcing_itp(δcalcite_file, SST_file, ∇T_file, Δt)
    
        δcalcite_itp, δcalcite_lims = δcalcite_interpolator(δcalcite_file);
        SST_itp, SST_lims = SST_interpolator(SST_file);
        ∇T_s_itp, ∇T_w_itp, ∇T_lims = ∇T_interpolator(∇T_file);
    
        age = minimum([δcalcite_lims[2],SST_lims[2],∇T_lims[2]]):-Δt:maximum([δcalcite_lims[1],SST_lims[1],∇T_lims[1]]);
    
        return age, δcalcite_itp, SST_itp, ∇T_s_itp, ∇T_w_itp

    end

    function δcalcite_interpolator(file)
        X = readdlm(file, ',');
        #offset = sum(X[1:4,2])./4.0
        itp = LinearInterpolation(X[:,1], X[:,2]);
        age_lim = [minimum(X[:,1]),maximum(X[:,1])]
        return itp, age_lim
    end

    function SST_interpolator(file)
        X = readdlm(file, ',');
        itp = LinearInterpolation(X[:,1], X[:,2]);
        age_lim = [minimum(X[:,1]),maximum(X[:,1])]
        return itp, age_lim
    end

    function ∇T_interpolator(file)
        X = readdlm(file, ',');
        itp1 = LinearInterpolation(X[:,1], X[:,2]);
        itp2 = LinearInterpolation(X[:,1], X[:,3]);
        age_lim = [minimum(X[:,1]),maximum(X[:,1])]
        return itp1, itp2, age_lim
    end