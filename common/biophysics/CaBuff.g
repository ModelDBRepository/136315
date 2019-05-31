//==================================================================
//                Simple Calcium Concentration Handling 
//==================================================================

float tau_CaClearance = 0.001
float shell_thick = 20e-9
float B_Ca_GP_conc = 5.2e-12

// Concentration object will keep track of I(Ca2+) and apply buffering.
function make_Ca_GP_conc
    if (({exists Ca_GP_conc}))
        return
    end
    create Ca_concen Ca_GP_conc
    setfield Ca_GP_conc    \
        tau        {tau_CaClearance}   \
        B        {B_Ca_GP_conc}         \
        Ca_base 5e-05          //Units in mM, so = 50 nM.
end

// Nernst object keeps track of Calcium reversal potential
function make_Ca_GP_nernst
    if (({exists Ca_GP_nernst}))
        return
    end
    create nernst Ca_GP_nernst
    setfield Ca_GP_nernst    \
        Cout    2    \    //external Ca2+ conc
        Cin    5e-5    \    //baseline internal Ca2+ conc 50 nM
        T    32    \    //temp in Celsius
        valency    2    \    //divalent
        scale     1        //E in volts
end
