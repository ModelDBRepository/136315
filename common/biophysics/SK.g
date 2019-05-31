//==================================================================
//        SK channel from Volker Steuber's DCN neuron model, 
//          Steuber V, Schultheiss NW, Silver RA, De Schutter E, Jaeger D
//            J Comput Neurosci 2010.
//    
//        Steuber version modified to reflect Hill fits in the following:
//        Hirschberg et al, 1999: Biophys J. 77: 1905-13. 
//        Keen et al, 1999: J. Neurosci 19: 8830-38.
//        Tau-Ca equation made by Volker based on Hirschberg et al, 1998:
//            JGP 111: 565-581.
//==================================================================    

// --> kinetics for room temperature

float zpower_SK       = 1    
float EC50_SK         = 0.00035   // SI unit = mM; default = 350 nM.
float hillslope_SK    = 4.6       // Hirschberg et al, 1999
float taumin_SK       = 0.008     // fastest tau-act in saturating Ca2+
float taumax_SK       = 0.076     // deactivation time constant in 0 Ca2+
float CaSat_SK        = 0.005     // calcium conc at which tau-act reaches max.

float dq10_SK = 2
function make_SK_GP
    if ({exists K_ahp_GP})
        return
    end
 
    float cxmin, cxmax, cxdivs, cdx
    float taum, minf
    float hillslope = {hillslope_SK}    // Hirschberg et al, 1999
    float taumax = {taumax_SK}          // deactivation tau in 0 Ca2+
    float taumin = {taumin_SK}          // max rate in saturating Ca2+
    float caSat = {CaSat_SK}            // calcium conc at which tauact hits max
    float tauK = ({taumax} - {taumin}) / {caSat}
    // NOTE: genesis SI concentration units = mols / m^3 = millimolar!
    create tabchannel K_ahp_GP
    setfield K_ahp_GP Ek {EK} Gbar {G_K_ahp_GP}  \
        Xpower 0 Ypower 0 Zpower {zpower_SK} 
    cxmin = 0.00001    // 10 nM
    cxmax = 0.06       // 60 microM
                       // Conc-dependence grid will have 6000 points spanning 
                       //    1 nM to 6 microM with resolution of 1 nM.
    cxdivs = 5999      // Have to use high-resolution for channel setup because
                       // most of G-Ca curve falls within the first 1 microM.
    cdx = (cxmax - cxmin)/cxdivs
    call K_ahp_GP TABCREATE Z {cxdivs} {cxmin} {cxmax}
    x = cxmin

    for (i = 0; i <= {cxdivs}; i = i + 1)
    
        if (x < {caSat})
              taum = ({taumax} - x*{tauK})/{dq10_SK}
        else
              taum = {taumin}/{dq10_SK}
        end
        minf = {pow {x} {hillslope}} / ({pow {x} {hillslope}} + {pow {EC50_SK} {hillslope}})

        setfield K_ahp_GP Z_A->table[{i}] {taum}
        setfield K_ahp_GP Z_B->table[{i}] {minf}

        x = x + cdx
      end
    
    tweaktau K_ahp_GP Z
    call K_ahp_GP TABFILL Z 6000 0
    setfield K_ahp_GP Z_A->calc_mode {NO_INTERP}
    setfield K_ahp_GP Z_B->calc_mode {NO_INTERP}
end
