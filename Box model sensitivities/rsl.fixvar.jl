#module fixvar #functions for fixed variables

    function fixvar() # get fixed parameters
        
        P_p = 1.25E12; #Precipitation rate (units: m^3 yr^-1)
        SA_p = 36.2; #Gibraltar Strait inflow salinity
        δA_p = 1.0; #Gibraltar Strait inflow oxygen isotope ratio
        X_p = 1.3E12; #Excess of evaporation over freshwater input (units: m^3 yr^-1)
        R_p = 0.45E12; #Runoff rate (units: m^3 yr^-1)
        A_p = 23.0E12; #Atlantic inflow through Strait of Gibraltar (units: m^3 yr^-1)
        B_p = 0.2E12; #Black Sea input (units: m^3 yr^-1)
        M_p = 21312*60*60*24*365; #Monsoon input (units: m^3)
        δB_p = -9.0; #mean d18O for B
        δR_p = -6.0; #mean d18O for R during non-monsoon times
        sinc_p = 1.0;
        δM_p = -10.0; #mean d18O for monsoon input

        return A_p, SA_p, δA_p, P_p, X_p, R_p, δR_p, B_p, δB_p, M_p, δM_p, sinc_p
    end

#end