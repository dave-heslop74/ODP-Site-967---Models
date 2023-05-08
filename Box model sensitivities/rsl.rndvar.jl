    function rand_r() #draw relative humidity
        r = 0.7.+0.05.*randn(Ntrials)./3.0
        return r
    end

    function rand_V() #draw wind speed 
        V = 7.5.+sl./120.0.+randn(Ntrials)./3.0
        return V
    end

    function rand_zsml() #draw depth of summer mixed layer
        #zsml = 30.0.+5.0.*randn(Ntrials)./3.0
        zsml = 50.0.+0.0.*randn(Ntrials)./3.0
        return zsml
    end

    function rand_δ() #draw δ values for R, P and B
        δR_w = δR_p.-1.0.+randn(Ntrials)./3.0 #winter runoff
        δR_s = δR_p.+1.0.+randn(Ntrials)./3.0 #summer runoff
        δP_w = copy(δR_w) #winter precipitation
        δP_s = copy(δR_s) #summer precipitation
        δB_w = δB_p.+randn(Ntrials)./3.0 #winter Black sea
        δB_s = δB_p.+randn(Ntrials)./3.0 #summer Black sea
        δM_m = δM_p.+randn(Ntrials)./3.0 #Monsoon runoff

        return δR_w, δR_s, δP_w, δP_s, δB_w, δB_s, δM_m
    end

    function rand_B() #draw Black Sea water flux (m^3 yr^-1)
        
        B = zeros(Ntrials)
        
        if sl <= 80
            B = (1.0.+0.1.*randn(Ntrials)./3.0).*B_p
        end
        
        return B
    end

    function rand_R() #draw Runoff water flux (m^3 yr^-1)
        R = (1.0.+0.1.*randn(Ntrials)./3.0).*R_p
        return R
    end

    function rand_T() #draw temperature parameters
        #println(length(T_s))
        
        if null_T
            T_s = T.+4.0 .- (∇T_s .+randn(Ntrials)./3.0).*(sl./120.0) #summer temperature
            #T_w = T.-4.0 .- (∇T_w .+randn(Ntrials)./3.0).*(sl./120.0) #winter temperature
            T_w = T.-2.0 .- (∇T_w .+randn(Ntrials)./3.0).*(sl./120.0) #winter temperature (Amies 2019)
        else
            T_s = μs_itp(age_i).+randn(Ntrials).*σs_itp(age_i); #Summer temperature based on estimated sensitivity
            T_w = μw_itp(age_i).+randn(Ntrials).*σw_itp(age_i); #Winter temperature based on estimated sensitivity
        end

        Tc_w = T_w.-(5.0.+randn(Ntrials)./3.0);
        Tc_s = T_s.-(5.0.+randn(Ntrials)./3.0);
        
        ΔT_s = 0.5.+(1.0.+randn(Ntrials)./3.0).*(sl./120.0) #Summer air-sea temperature difference
        ΔT_w = 1.5.+randn(Ntrials)./3.0 #Winter air-sea temperature difference
        
        Ta_s = T_s.-ΔT_s #Summer atmospheric temperature at 10 m above sea surface
        Ta_w = T_w.-ΔT_w #Winter atmospheric temperature at 10 m above sea surface

        Tssth = copy(T_w) #temperature of summer sub-thermocline
        Tint = T_w.-(1.0.+randn(Ntrials)./3.0) #temperature of intermediate waters
        
        TA = Tint.+randn(Ntrials)./3.0 #temperature of Atlantic inflow

        return T_s, T_w, Tc_s, Tc_w, Ta_s, Ta_w, Tssth, Tint, TA
    end 