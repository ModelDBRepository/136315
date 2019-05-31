//==================================================================
//                       HVA Ca2+ Channels 
//    Voltage-dependent activation from GP data:  
//        Surmeier Seno and Kitai 1994
//             J Neurophysio. 71(3):1272-1280
//==================================================================

// --> kinetics for 32 degrees C (Q10=2.5 assumption)

float npower_CaHVA     = 1
float Vhalfn_CaHVA     = -0.02
float Kn_CaHVA         = 0.007 
float taun_CaHVA    = 0.0002

float dq10_CaHVA    = 1
    
function make_Ca_HVA_GP
    if (({exists Ca_HVA_GP}))
        return
    end
    int ndivs, i
    float x, y
    create tabchannel Ca_HVA_GP
    setfield Ca_HVA_GP Ek {ECa} Gbar {G_Ca_HVA_GP} Ik 0 Gk 0 \
        Xpower {npower_CaHVA} Ypower 0 Zpower 0

    //first setup voltage-dependent activation
    float tau =  {taun_CaHVA} / {dq10_CaHVA}
    float K = -1*{Kn_CaHVA}
    float V0 = {Vhalfn_CaHVA}
    setuptau Ca_HVA_GP X \
        {tau} {tau*1e-6} 0 0 1e6 \
        1 0 1 {-1.0*V0} {K} -range {xmin} {xmax} 
    call Ca_HVA_GP TABFILL X 6000 0
    setfield Ca_HVA_GP X_A->calc_mode {NO_INTERP}
    setfield Ca_HVA_GP X_B->calc_mode {NO_INTERP}
end
