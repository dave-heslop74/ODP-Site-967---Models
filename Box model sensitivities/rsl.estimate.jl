
    function T_to_qsat(T,p) # convert temperature into saturation mixing ratio
        TK = T.+273.15 #temperature in Kelvin
        esat = 1E2.*exp.(55.17 .- 6803.0./TK - 5.07.*log.(TK))
        qsat = 18.0153./28.965 * esat./(p.-esat)

        return qsat
    end

    function T_to_fluxes(SST,Ta) #convert T into water mass fluxes
    
        Ntrials = length(V)
        L = (2500.84.-2.34.*SST).*1.0E3 #Latent heat of vaporization
    
        qsat_SST = T_to_qsat(SST,p)
        qsat_Ta = T_to_qsat(Ta,p)
        Z = qsat_Ta ./ qsat_SST
    
        QL = ρa.*L.*C.*V.*(qsat_SST-r.*qsat_Ta)
        E = QL.*AM.*1.26E-2 #Evaporation
        E = E.*Δt #season only considers portion of year

        P = 0.4.*E #Precipitation
        R = (1.0.+0.1.*randn(Ntrials)./3.0).*R_p.*P./P_p
        M = M_p .* ones(Ntrials) .* Mfactor
        
        if season == 's' || season == 'm'
            Bnew = B_p.*Δt.+sinc_p.*(B.-B_p)
        else
            Bnew = B_p.*Δt.+(1.0.-sinc_p).*(B.-B_p)
        end
        
        X = E.-P.-R.-Bnew.-M
        
        return E, P, R, Bnew, M, X
    end

     function sl_to_inflow(SA_p,δA_p,sl) #Altantic inflow properties as function of sea-level
        SA = SA_p .+ sl./120.0
        δA = δA_p .+ sl.*0.009 
        
        return SA, δA 
    end

     function T_to_δEequ(T) 
        αevap = exp.(1.137./(T.+273.15).^2.0.*1E3.-0.4156./(T.+273.15).-2.0667E-3)
        δEequ = -log.(αevap).*1E3

        return δEequ, αevap
    end

     function T_to_δatm(T,δR) 
        α = exp.(1.137./(T.+273.15).^2.0.*1E3.-0.4156./(T.+273.15).-2.0667E-3)
        δatm = δR.-log.(α).*1E3
        
        return δatm
    end

     function δ_to_δcalcite(δM,T)
        αcalc = exp.(2.78.*(T.+273.15).^-2 .* 1.0E3 .- 3.39E-3)
        δcalcite = 1.0E3 .* log.(αcalc) .+ (δM.-30.92)./1.03092

        return δcalcite
    end

     function δM_to_δE(SST,Ta,δM,αevap,r,p,δatm)
        
        Z = T_to_qsat(Ta,p)./T_to_qsat(SST,p)
        δM0 = δM .* 1.0E-3
        δatm0 = δatm .* 1.0E-3
        δE0 = (1.0.+δM0)./ αevap .- Z .* r .* (1 .+ δatm0)
        δE0 ./= (1.0 .- Z .* r).*1.0142
        δE = (δE0.-1) .* 1E3
    
        return δE
    end

     function atlantic_flux()
        
        #SA_p = present day salinity of Gibralter inflow
        #δA_p = present day oxygen isotope composition of Gibralter inflow
        #A_p = present day inflow from the Atlantic
        #X_p = present day excess evaporation
        #sl = sea-level
        #X_s = summer excess of evaporation
        #X_m = monsoon excess of evaporation
        #X_w = winter excess of evaporation
        #zsml = depth of summer mixed layer
        #zmodel = depth of entire model
        #AM = Area of mediterranean
    
        SA = SA_p.+sl./120
        δA = δA_p.+sl.*0.009
        Φ = 1.0.-0.5.*sl./120

        X = X_s.+X_m.+X_w
        γ = X./X_p
        γ[γ.<1] .= 1;
        Ω = γ.^(1.0./3.0)
    
        A = Ω .* Φ .* A_p #Atlantic inflow flux
    
        α = (A_p./X_p.*Φ.*Ω-γ)./(A_p./X_p.-1.0).*(SA_p./SA)
        zssth = α.*zmodel.-zsml

        return A, SA, δA, zssth
    end   

    function estimate_summer(Δt)

        #T = SST of Mediterranean
        #Ta = atm temperature above Mediterranean
        #δM = δ of Mediterranean
        #SM = Salinity of Mediterranean
        #δA = δ of Atlantic
        #SA = Salinity of Atlantic
        #A = Altantic flux
        #X = Excess evaporation flux
        #R = Freshwater runoff flux
        #δR = δ of freshwater runoff
        #B = Black Sea input flux
        #δB = δ of Black Sea input
        #P = Precipitation flux
        #δP = δ of precipitation
        #zsml = thickness of sml 
        #zssth = thickness of ssth        
        #AM = area of Mediterranean
        #Δt_s = summer time step (fraction of year)
        #Δt_w = winter time step (fraction of year)
        
        #initial estimate of δE
        δEequ, αevap = T_to_δEequ(T_s)
        δE_s = δwml .+ δEequ

        #Estimate SSTH based on Atlantic / Mediterranean mixture
        VA = A.*Δt.*zssth./(zsml.+zssth) #Volume of Atlantic water in Box
        VM = zssth.*AM.-VA #Volume of water in Mediterranean box
        δssth = (δA.*VA .+ δwml.*VM)./(VA.+VM) # mixture of δ
        Sssth = (SA.*VA .+ Swml.*VM)./(VA.+VM) # mixture of salinities

        #Initial SML will be the same as SSTH
        δsml = copy(δssth) # mixture of δ
        δsml0 = copy(δssth) # mixture of δ for iterative estimate
        Ssml = copy(Sssth) # mixture of salinities

        #Add additional δ-surface fluxes to SML 
        δsml .*= zsml.*AM
        δsml .+= (R_s.*δR_s + B_s.*δB_s + P_s.*δP_s - E_s.*δE_s)
        δsml ./= (zsml.*AM.-X_s)

        #Add additional Salinity-surface fluxes to SML 
        Ssml .*= zsml.*AM
        Ssml ./= (zsml.*AM.-X_s)

        Tc = T_s.-(5.0.+randn(length(T_s))./3.0)
        δatm = T_to_δatm(Tc,δR_s)

        for i in 1:3 #Iterate δsml with expression for δE to find equilibrium solution
            δE_s = δM_to_δE(T_s,Ta_s,δsml,αevap,r,p,δatm)
        
            #Estimate δ for sml including additional fluxes
            δsml = δsml0 .* zsml .* AM
            δsml .+= (R_s.*δR_s .+ B_s.*δB_s .+ P_s.*δP_s .- E_s.*δE_s)
            δsml ./= (zsml.*AM .- X_s)
        end
    
        return δsml, Ssml, δssth, Sssth

    end
   

     function estimate_monsoon(Δt)

        #T = SST of Mediterranean
        #Ta = atm temperature above Mediterranean
        #δM = δ of Mediterranean
        #SM = Salinity of Mediterranean
        #δA = δ of Atlantic
        #SA = Salinity of Atlantic
        #A = Altantic flux
        #X = Excess evaporation flux
        #R = Freshwater runoff flux
        #δR = δ of freshwater runoff
        #B = Black Sea input flux
        #δB = δ of Black Sea input
        #P = Precipitation flux
        #δP = δ of precipitation
        #zsml = thickness of sml 
        #zssth = thickness of ssth        
        #AM = area of Mediterranean
        #Δt_s = summer time step (fraction of year)
        #Δt_w = winter time step (fraction of year)
        
        #estimate thickness of all freshwater contributions
        zmll = zsml .- zmul #set thickness of lower monsoon layer

        #initial estimate of δE
        δEequ, αevap = T_to_δEequ(T_m)
        δE_m = δsml + δEequ

        #Estimate MSTH based on Atlantic / Mediterranean mixture
        VA = A.*Δt.*zssth./(zsml.+zssth) #Volume of Atlantic water in Box
        VM = zssth.*AM.-VA #Volume of water in Mediterranean box
        δmsth = (δA.*VA .+ δssth.*VM)./(VA.+VM) # mixture of δ
        Smsth = (SA.*VA .+ Sssth.*VM)./(VA.+VM) # mixture of salinities

        #Estimate MLL based on Atlantic / Mediterranean mixture
        VA = A.*Δt.*zmll./(zsml.+zssth) #Volume of Atlantic water in Box
        VM = zmll.*AM.-VA #Volume of water in Mediterranean box
        δmll = (δA.*VA .+ δssth.*VM)./(VA.+VM) # mixture of δ
        Smll = (SA.*VA .+ Sssth.*VM)./(VA.+VM) # mixture of salinities

        #estimate MUL based on Atlantic / Mediterranean mixture
        VA = A.*Δt.*zmul./(zsml.+zssth) #Volume of Atlantic water in Box
        VM = zmul.*AM.-X_s.-VA #Volume of water in Mediterranean box
        δmul = (δA.*VA .+ δsml.*VM)./(VA.+VM) # mixture of δ
        δmul0 = copy(δmul) # mixture of δ for iterative estimate
        Smul = (SA.*VA .+ Ssml.*VM)./(VA.+VM) # mixture of salinities

        #Add additional δ-surface fluxes to MUL 
        δmul .*= zmul.*AM.-X_s
        δmul .+= (R_m.*δR_m .+ B_m.*δB_m .+ P_m.*δP_m .+ M_m.*δM_m .- E_m.*δE_m)
        δmul ./= (zmul.*AM.-X_s.-X_m)

        #Add additional Salinity-surface fluxes to MUL 
        Smul .*= zmul.*AM.-X_s
        Smul ./= zmul.*AM.-X_s.-X_m

        Tc = T_m.-(5.0.+randn(length(T_m))./3.0)
        δatm = T_to_δatm(Tc,δR_m)

        for i in 1:3 #Iterate δsml with expression for δE to find equilibrium solution
            δE_m = δM_to_δE(T_m,Ta_m,δmul,αevap,r,p,δatm)
        
            #Estimate δ for sml including additional fluxes
            δmul = δmul0 .* (zmul .* AM .- X_s)
            δmul .+= (R_m.*δR_m .+ B_m.*δB_m .+ P_m.*δP_m .+ M_m.*δM_m .- E_m.*δE_m)
            δmul ./= (zmul.*AM .- X_s .- X_m)
        end
    
        return δmul, Smul, δmll, Smll, δmsth, Smsth, zmll

    end


    function estimate_summer2winter(Δt)

        #T = Summer SST of Mediterranean
        #Ta = Summer atm temperature above Mediterranean
        #δM = δ of Mediterranean
        #SM = Salinity of Mediterranean
        #δA = δ of Atlantic
        #SA = Salinity of Atlantic
        #A = Altantic flux
        #X = Excess evaporation flux
        #R = Freshwater runoff flux
        #δR = δ of freshwater runoff
        #B = Black Sea input flux
        #δB = δ of Black Sea input
        #P = Precipitation flux
        #δP = δ of precipitation
        #zsml = thickness of sml 
        #zssth = thickness of ssth        
        #AM = area of Mediterranean
        #Δt_s = summer time step (fraction of year)
        #Δt_w = winter time step (fraction of year)
        
        #initial estimate of δE for iterative process
        # Previous year was not monsoon, so need to convert summer boxes into fake monsoon boxes
       

        δEequ, αevap = T_to_δEequ(T_w)
        δE_w = (zsml.*δsml.+zssth.*δmsth)./(zsml.+zssth).+δEequ

        #Initial mixing of SML & SSTH to form WML
        δwml = δsml.*zsml.*AM.-X_s
        δwml .+= (δmsth.*zssth.*AM)
        δwml ./= (zsml.*AM.+zssth.*AM.-X_s)
        
        Swml = Ssml.*(zsml.*AM.-X_s)
        Swml .+= (Sssth.*zssth.*AM)
        Swml ./= (zsml.*AM.+zssth.*AM.-X_s)

        #Estimate MSML based on Atlantic / Mediterranean mixture
        VA = A.*Δt #Volume of Atlantic water in Box
        VM = (zsml.+zssth).*AM.-X_s.-VA #Volume of water in Mediterranean box
        δwml = (δA.*VA .+ δwml.*VM)./(VA.+VM) # mixture of δ
        δwml0 = copy(δwml) # mixture of δ for iterative estimate
        Swml = (SA.*VA .+ Swml.*VM)./(VA.+VM) # mixture of salinities

        #Add additional Salinity-surface fluxes to WML 
        Swml .*= (zsml.+zssth).*AM.-X_s
        Swml ./= (zsml.+zssth).*AM.-X_s.-X_w

        #Add additional δ-surface fluxes to WML
        δwml .*= (zsml.+zssth).*AM.-X_s
        δwml .+= (R_w.*δR_w + B_w.*δB_w + P_w.*δP_w .- E_w.*δE_w)
        δwml ./= (zsml.+zssth).*AM.-X_s.-X_w

        Tc = T.-(5.0.+randn(length(T))./3.0)
        δatm = T_to_δatm(Tc,δR_w)

        for i in 1:3 #Iterate δwml with expression for δE to find equilibrium solution
            δE_w = δM_to_δE(T_w,Ta_w,δwml,αevap,r,p,δatm)
        
            δwml = δwml0 .* ((zsml.+zssth) .* AM .- X_s)
            δwml .+= (R_w.*δR_w .+ B_w.*δB_w .+ P_w.*δP_w .- E_w.*δE_w)
            δwml ./= ((zsml.+zssth) .* AM .- X_s .- X_w)
        end
    
        return δwml, Swml

    end


    function estimate_winter(Δt)

        #T = Summer SST of Mediterranean
        #Ta = Summer atm temperature above Mediterranean
        #δM = δ of Mediterranean
        #SM = Salinity of Mediterranean
        #δA = δ of Atlantic
        #SA = Salinity of Atlantic
        #A = Altantic flux
        #X = Excess evaporation flux
        #R = Freshwater runoff flux
        #δR = δ of freshwater runoff
        #B = Black Sea input flux
        #δB = δ of Black Sea input
        #P = Precipitation flux
        #δP = δ of precipitation
        #zsml = thickness of sml 
        #zssth = thickness of ssth        
        #AM = area of Mediterranean
        #Δt_s = summer time step (fraction of year)
        #Δt_w = winter time step (fraction of year)
        
        #initial estimate of δE for iterative process
        # Previous year was not monsoon, so need to convert summer boxes into fake monsoon boxes
       

        δEequ, αevap = T_to_δEequ(T_w)
        δE_w = (zmul.*δmul.+zmll.*δmll.+zssth.*δmsth)./(zmul.+zmll.+zssth).+δEequ

        #Initial mixing of MML & MSTH to form WML
        δwml = δmul.*zmul.*AM.-X_s.-X_m
        δwml .+= (δmll.*zmll.*AM)
        δwml .+= (δmsth.*zssth.*AM)
        δwml ./= (zmul.*AM.+zmll.*AM.+zssth.*AM.-X_s.-X_m)
        
        Swml = Smul.*(zmul.*AM.-X_s.-X_m)
        Swml .+= (Smll.*zmll.*AM)
        Swml .+= (Smsth.*zssth.*AM)
        Swml ./= (zmul.*AM.+zmll.*AM.+zssth.*AM.-X_s.-X_m)

        #Estimate MSML based on Atlantic / Mediterranean mixture
        VA = A.*Δt #Volume of Atlantic water in Box
        VM = (zmul.+zmll.+zssth).*AM.-X_s.-X_m.-VA #Volume of water in Mediterranean box
        δwml = (δA.*VA .+ δwml.*VM)./(VA.+VM) # mixture of δ
        δwml0 = copy(δwml) # mixture of δ for iterative estimate
        Swml = (SA.*VA .+ Swml.*VM)./(VA.+VM) # mixture of salinities

        #Add additional Salinity-surface fluxes to WML 
        Swml .*= (zmul.+zmll.+zssth).*AM.-X_s.-X_m
        Swml ./= (zmul.+zmll.+zssth).*AM.-X_s.-X_m.-X_w

        #Add additional δ-surface fluxes to WML
        δwml .*= (zmul.+zmll.+zssth).*AM.-X_s.-X_m
        δwml .+= (R_w.*δR_w + B_w.*δB_w + P_w.*δP_w .- E_w.*δE_w)
        δwml ./= (zmul.+zmll.+zssth).*AM.-X_s.-X_m.-X_w

        Tc = T.-(5.0.+randn(length(T))./3.0)
        δatm = T_to_δatm(Tc,δR_w)

        for i in 1:3 #Iterate δwml with expression for δE to find equilibrium solution
            δE_w = δM_to_δE(T_w,Ta_w,δwml,αevap,r,p,δatm)
        
            δwml = δwml0 .* ((zmul.+zmll.+zssth) .* AM .- X_s .- X_m)
            δwml .+= (R_w.*δR_w .+ B_w.*δB_w .+ P_w.*δP_w .- E_w.*δE_w)
            δwml ./= ((zmul.+zmll.+zssth) .* AM .- X_s .- X_m .- X_w)
        end
    
        return δwml, Swml

    end

    function estimate_years(Ntrials,Nyrs,annual_forcing)

        #predefine random forcing
        global V = rand_V(); #wind speed
        global r = rand_r(); #relative humidity
        global zsml = rand_zsml(); #depth of summer mixed layer
        global B = rand_B(); #Black sea flux
        global R = rand_R(); #Runoff flux
        global δR_w, δR_s, δP_w, δP_s, δB_w, δB_s, δM_m = rand_δ(); #seasonal freshwater δ values
        global T_s, T_w, Tc_s, Tc_w, Ta_s, Ta_w, Tssth, Tint, TA = rand_T(); #Temperature forcing values
        global T_m, Ta_m, δR_m, δP_m, δB_m = copy(T_s).+ΔT_m,copy(Ta_s).+ΔT_m,copy(δR_s), copy(δP_s), copy(δB_s); #Assign monsoon parameters

        #predefine estimated SUMMER fluxes
        global season, Δt, Mfactor = 's',Δt_s,Mf_s #Summer parameters
        global E_s, P_s, R_s, B_s, M_s, X_s = T_to_fluxes(T_s,Ta_s); #Summer water fluxes
        
        #predefine estimated MONSOON fluxes
        season, Δt, Mfactor = 'm',Δt_m,Mf_m #Monsoon parameters
        global E_m, P_m, R_m, B_m, M_m, X_m = T_to_fluxes(T_m,Ta_m); #Monsoon water fluxes    
        if Δt==0
            E_m, P_m, R_m, B_m, M_m, X_m = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        end

        #predefine estimated WINTER fluxes
        season, Δt, Mfactor = 'w',Δt_w,Mf_w #Winter parameters
        global E_w, P_w, R_w, B_w, M_s, X_w = T_to_fluxes(T_w,Ta_w); #Winter water fluxes

        #predefine estimated ATLANTIC flux
        global A, SA, δA, zssth = atlantic_flux();

        if spinup
            global δwml = copy(δA.*ones(Ntrials)); #initial δwml value for iterative estimation
            global Swml = copy(SA.*ones(Ntrials)); #initial Swml value for iterative estimation
        end

        for i in 1:Nyrs
            
            if annual_forcing #sample forcing each year
                V = rand_V(); #wind speed
                r = rand_r(); #relative humidity
                zsml = rand_zsml(); #depth of summer mixed layer
                B = rand_B(); #Black sea flux
                R = rand_R(); #Runoff flux
                δR_w, δR_s, δP_w, δP_s, δB_w, δB_s, δM_m = rand_δ(); #seasonal freshwater δ values
                T_s, T_w, Tc_s, Tc_w, Ta_s, Ta_w, Tssth, Tint, TA = rand_T(); #Temperature forcing values
                T_m, Ta_m, δR_m, δP_m, δB_m = copy(T_s).+ΔT_m,copy(Ta_s).+ΔT_m,copy(δR_s), copy(δP_s), copy(δB_s); #Assign monsoon parameters
            
                #redefine estimated SUMMER fluxes
                season, Δt, Mfactor = 's', Δt_s, Mf_s #Summer parameters
                E_s, P_s, R_s, B_s, M_s, X_s = T_to_fluxes(T_s,Ta_s); #Summer water fluxes
                
                #redefine estimated MONSOON fluxes
                season, Δt, Mfactor = 'm', Δt_m, Mf_m #Monsoon parameters
                if Δt>0 
                    E_m, P_m, R_m, B_m, M_m, X_m = T_to_fluxes(T_m,Ta_m); #Monsoon water fluxes
                else
                    E_m, P_m, R_m, B_m, M_m, X_m = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                end    

                #redefine estimated WINTER fluxes
                season, Δt, Mfactor = 'w',Δt_w, Mf_w #Winter parameters
                E_w, P_w, R_w, B_w, M_s, X_w = T_to_fluxes(T_w,Ta_w); #Winter water fluxes

                #redefine estimated ATLANTIC flux
                A, SA, δA, zssth = atlantic_flux();
            end   
            
            global δsml, Ssml, δssth, Sssth = estimate_summer(Δt_s);
            global δmul, Smul, δmll, Smll, δmsth, Smsth, zmll = estimate_monsoon(Δt_m);

            if Δt_m>0
                δwml, Swml = estimate_winter(Δt_w);
            else
                δwml, Swml = estimate_summer2winter(Δt_w);
            end
        end
        
        δcalcite_w = δ_to_δcalcite(δwml,T_w)
        δcalcite_s = δ_to_δcalcite(δsml,T_s)
        δcalcite_m = δ_to_δcalcite(δmul,T_m)
        
        δcalcite_m0 = δcalcite_m.-(δcalcite_m.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999
        T_m0 = T_m.-(T_m.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999
        δmul0 = δmul.-(δmul.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999
        δmll0 = δmll.-(δmll.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999        
        δmsth0 = δmsth.-(δmsth.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999         
        Smul0 = Smul.-(Smul.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999
        Smll0 = Smll.-(Smll.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999        
        Smsth0 = Smsth.-(Smsth.*(1-ceil(Δt_m))).-(1-ceil(Δt_m)).*9999         
        
        return δsml, Ssml, δssth, Sssth, δmul0, Smul0, δmll0, Smll0, δmsth0, Smsth0, δwml, Swml, T_w, T_s, T_m0, δcalcite_w, δcalcite_s, δcalcite_m0
    
    end