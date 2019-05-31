//==================================================================
//              Kdr Kv2 
//              (Kv2.1) slow activating
//              Created based on GP data:
//                Baranuskas, Tkatch and Surmeier
//                  1999, J Neurosci 19(15):6394-6404
//==================================================================

// --> kinetics for 32 degrees C:
float npower_Kv2    = 4
float Vhalfn_Kv2    = -0.018
float Kn_Kv2        = 0.0091
float taunmin_Kv2   = 0.0001
float taunmax_Kv2   = 0.03
float Ktaun1_Kv2    = 0.02174
float Ktaun2_Kv2    = -0.01391

float hpower_Kv2    = 1
float hmin_Kv2      = 0.2
float V0h_Kv2       = -0.02
float Kh_Kv2        = -0.01
float tauhmin_Kv2   = 3.4
float tauhmax_Kv2   = 3.4
float V0tauh_Kv2    = 0        // irrelevant while tauhmin == tauhmax
float Ktauh1_Kv2    = 0.01    // irrelevant while tauhmin == tauhmax
float Ktauh2_Kv2    = -0.01 // irrelevant while tauhmin == tauhmax

float dq10_Kv2      = 1

function make_Kv2_GP
    if (({exists Kv2_GP}))
        return
    end
    create tabchannel Kv2_GP
    setfield Kv2_GP Ek {EK} Gbar {G_Kv2_GP} Ik 0 Gk 0\
        Xpower {npower_Kv2} Ypower {hpower_Kv2} Zpower 0

    float Vhalfn  = {Vhalfn_Kv2}    // True Vhalf for channel activation
    float Kn      = {Kn_Kv2}
    float taunmax = {taunmax_Kv2} / {dq10_Kv2}
    float taunmin = {taunmin_Kv2} / {dq10_Kv2}
    float K1tau   = {Ktaun1_Kv2}
    float K2tau   = {Ktaun2_Kv2}

    float V0n, ninf, taun, alpha, beta
    V0n = {Vhalfn} + ({Kn} * {log {(1 / {pow 0.5 {1/{npower_Kv2}}}) - 1}})
    echo "Kv2 V0n: " {V0n}
    //V0n is Vhalf for each individual n gate.
    call Kv2_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        ninf = 1 / (1 + {exp { ({V0n} - x) / {Kn} } } )
        taun = {taunmin} + (({taunmax} - {taunmin}) / ({exp { ({V0n} - x) / {K1tau} } } + {exp {({V0n} - x) / {K2tau} }}))
        setfield Kv2_GP X_A->table[{i}] {taun}
        setfield Kv2_GP X_B->table[{i}] {ninf}
        x = x + dx
    end
    tweaktau Kv2_GP X
    call Kv2_GP TABFILL X 6000 0
    setfield Kv2_GP X_A->calc_mode {NO_INTERP}
    setfield Kv2_GP X_B->calc_mode {NO_INTERP}
        
    float V0h     = {V0h_Kv2}
    float Kh      = {Kh_Kv2}
    float hmin    = {hmin_Kv2}
    float tauhmax = {tauhmax_Kv2} / {dq10_Kv2}
    float tauhmin = {tauhmin_Kv2} / {dq10_Kv2}
    float Ktauh1  = {Ktauh1_Kv2}
    float Ktauh2  = {Ktauh2_Kv2}
    float V0tauh  = {V0tauh_Kv2} 
    float hinf, tauh
    call Kv2_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        hinf = {hmin} + ((1 - {hmin}) / (1 + {exp {({V0h} - x) / {Kh} }}))
        tauh = {tauhmin} + (({tauhmax} - {tauhmin}) / ({exp { ({V0tauh} - x) / {Ktauh1}}} + {exp {({V0tauh} - x) / {Ktauh2} }}))
        setfield Kv2_GP Y_A->table[{i}] {tauh}
        setfield Kv2_GP Y_B->table[{i}] {hinf}
        x = x + dx
    end
    tweaktau Kv2_GP Y
    call Kv2_GP TABFILL Y 6000 0
    setfield Kv2_GP Y_A->calc_mode {NO_INTERP}
    setfield Kv2_GP Y_B->calc_mode {NO_INTERP}
end
