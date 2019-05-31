// *************************************************************************
function make_GP_AMPA
    if (!({exists AMPA}))
        create synchan AMPA
    end
    setfield AMPA Ek {E_AMPA} tau1 {tauRise_AMPA} tau2 {tauFall_AMPA} \
        gmax {G_AMPA} frequency 0 
end

// *************************************************************************
function make_GP_GABA_striatum
    if (!({exists GABA}))
               create synchan GABA
    end
    setfield GABA Ek {E_GABA} tau1 {tauRise_GABA} tau2 {tauFall_GABA}  \
        gmax {G_GABA} frequency 0
end

// *************************************************************************
function make_GP_GABA_pallidum
    if (!({exists GABA_GP}))
               create synchan GABA_GP
    end
    setfield GABA_GP Ek {E_GABA} tau1 {tauRise_GABA_GP} \
        tau2 {tauFall_GABA_GP} gmax {G_GABA_GP} frequency 0
end

// *************************************************************************
function make_GP_syns
    make_GP_AMPA
    make_GP_GABA_striatum
    make_GP_GABA_pallidum    
end
