//==================================================================
//          Kdr Kv3
//          (Kv3.1/3.4 heteromultimers) fast activating,
//            incompletely inactivating.
//          From Surmeier's kv3sur.mod NEURON mechanism
//              written by Josh Held
//          Adapted to genesis by J.R. Edgerton, 2003
//          Modified to reflect new data in
//            Baranauskas et al (2003) by J.R. Edgerton, 2004.
//            Baranauskas, Tkatch and Surmeier
//                  1999, J Neurosci 19(15):6394-6404
//            Baranauskas, Tkatch, Nagata, Yeh & Surmeier 2003.
//                  Nat Neurosci 6: 258-66.
//==================================================================

// --> kinetics for 32 degrees C:

float npower_Kv3    = 4
float Vhalfn_Kv3    = -0.013    // Actual Vhalf
float Kn_Kv3        = 0.0078    // Yields K = 6 mV with Xpower = 4
float taunmin_Kv3   = 0.0001    // 32 degrees C
float taunmax_Kv3   = 0.014        // 32 degrees C
float Ktaun1_Kv3    = -0.012
float Ktaun2_Kv3    = -0.013

float hpower_Kv3    = 1
float hmin_Kv3      = 0.6
float V0h_Kv3       = -0.02
float Kh_Kv3        = -0.010
float tauhmin_Kv3   = 0.007    
float tauhmax_Kv3   = 0.033    
float V0tauh_Kv3    = 0
float Ktauh1_Kv3    = 0.01
float Ktauh2_Kv3    = -0.01

float dq10_Kv3      = 1
function make_Kv3_GP
    if (({exists Kv3_GP}))
        return
    end
    create tabchannel Kv3_GP
    setfield Kv3_GP Ek {EK} Gbar {G_Kv3_GP} Ik 0 Gk 0\
        Xpower {npower_Kv3} Ypower {hpower_Kv3} Zpower 0
    
    float Vhalfn = {Vhalfn_Kv3}    // True Vhalf for channel activation
    float Kn = {Kn_Kv3}
    float taunmax = {taunmax_Kv3} / {dq10_Kv3}
    float taunmin = {taunmin_Kv3} / {dq10_Kv3}
    float K1tau = {Ktaun1_Kv3}
    float K2tau = {Ktaun2_Kv3}

    float V0n, ninf, taun, alpha, beta
    V0n = {Vhalfn} + ({Kn} * {log {(1 / {pow 0.5 {1/{npower_Kv3}}}) - 1}})
    echo "Kv3 V0n: " {V0n}
    //V0n is Vhalf for each individual n gate.
    
    call Kv3_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        ninf = 1 / (1 + {exp { ({V0n} - x) / {Kn} } } )
        taun = {taunmin} + (({taunmax} - {taunmin}) / ({exp { ({V0n} - x) / {K1tau} } } + {exp {-({V0n} - x) / {K2tau} }}))
        setfield Kv3_GP X_A->table[{i}] {taun}
        setfield Kv3_GP X_B->table[{i}] {ninf}
        x = x + dx
    end
    tweaktau Kv3_GP X
    call Kv3_GP TABFILL X 6000 0
    setfield Kv3_GP X_A->calc_mode {NO_INTERP}
    setfield Kv3_GP X_B->calc_mode {NO_INTERP}
        
    float V0h     =     {V0h_Kv3}
    float Kh      =      {Kh_Kv3}
    float hmin    =    {hmin_Kv3}
    float tauhmax = {tauhmax_Kv3} / {dq10_Kv3}
    float tauhmin =    {tauhmin_Kv3} / {dq10_Kv3}
    float Ktauh1  = {Ktauh1_Kv3}
    float Ktauh2  = {Ktauh2_Kv3}
    float V0tauh  =    {V0tauh_Kv3} 
    float hinf, tauh
    call Kv3_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        hinf = {hmin} + ((1 - {hmin}) / (1 + {exp {({V0h} - x) / {Kh} }}))
        tauh = {tauhmin} + (({tauhmax} - {tauhmin}) / ({exp { ({V0tauh} - x) / {Ktauh1}}} + {exp {({V0tauh} - x) / {Ktauh2} }}))
        setfield Kv3_GP Y_A->table[{i}] {tauh}
        setfield Kv3_GP Y_B->table[{i}] {hinf}
        x = x + dx
    end

    tweaktau Kv3_GP Y
    call Kv3_GP TABFILL Y 6000 0
    setfield Kv3_GP Y_A->calc_mode {NO_INTERP}
    setfield Kv3_GP Y_B->calc_mode {NO_INTERP}
end
