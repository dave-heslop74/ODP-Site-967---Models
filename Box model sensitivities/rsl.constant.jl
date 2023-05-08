function constants() # get constant parameters
    C = 1.15E-3; #Exchange Coefficient (units: ???)
    ρa = 1.2; #Density of air (units: kg m^-3)
    AM = 2.5E12; #Area of Mediterranean (units m^2)
    p = 1.012E5; #Total pressure at mean dea-level (units: Pa)

    return C, ρa, AM, p
end