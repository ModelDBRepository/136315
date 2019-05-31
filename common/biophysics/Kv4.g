//==================================================================
//                KA Kv4 fast and KA Kv4 slow
//                  Low voltage-activated
//                  Created based on GP data:
//                  Tkatch, Baranauskas and Surmeier
//                    2000, J Neurosci 20(2):579-588
//                  Reflects a mix of Kv4.1, Kv4.2, Kv4.3
//                  Modified by J. R. Edgerton 02/2004
//==================================================================

// --> n gate (activation/deactivation) is the same for Kv4-fast and Kv4-slow
// --> kinetics for 32 degrees C:

float npower_Kv4    = 4
float V0n_Kv4       = -0.049    // Yields Vhalf = -27 mV when Xpower = 4
float Kn_Kv4        = 0.0125    // Yields K = 9.6 mV when Xpower = 4
float taunmin_Kv4   = 0.00025
float taunmax_Kv4   = 0.007    
float Ktaun1_Kv4    = 0.029
float Ktaun2_Kv4    = -0.029

float hpower_Kv4    = 1
float V0h_Kv4       = -0.083
float Kh_Kv4        = -0.01    
float Ktauh1_Kv4    = 0.010
float Ktauh2_Kv4    = -0.010

// Only the inactivation time constants differ between Kv4f and Kv4s
float tauhmin_Kv4f  = 0.007
float tauhmax_Kv4f  = 0.021 
float tauhmin_Kv4s  = 0.050
float tauhmax_Kv4s  = 0.121

float dq10_Kv4      = 1

function make_Kv4_fast_GP
    if (({exists Kv4_fast_GP}))
        return
    end
    create tabchannel Kv4_fast_GP
    setfield Kv4_fast_GP Ek {EK} Gbar {G_Kv4_fast_GP} Ik 0 Gk 0\
        Xpower {npower_Kv4} Ypower {hpower_Kv4} Zpower 0
    float Kn        = {Kn_Kv4}
    float V0n       = {V0n_Kv4}
    float taunmax   = {taunmax_Kv4} / {dq10_Kv4}
    float taunmin   = {taunmin_Kv4} / {dq10_Kv4}
    float Ktaun1    = {Ktaun1_Kv4}
    float Ktaun2    = {Ktaun2_Kv4}
    float ninf, taun
    float Vhalfn = {V0n} - ({Kn} * {log {(1 / {pow 0.5 {1/{npower_Kv4}}}) - 1}})
    echo "Kv4f actual Vhalf: " {Vhalfn}
     call Kv4_fast_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        ninf = 1 / (1 + {exp {({V0n} - x) / {Kn} }})
        taun = taunmin + ({taunmax} - {taunmin}) / ({exp {({V0n} - x) / {Ktaun1}}} + {exp {({V0n} - x) / {Ktaun2}}}) 
        setfield Kv4_fast_GP X_A->table[{i}] {taun}
        setfield Kv4_fast_GP X_B->table[{i}] {ninf}
        x = x + dx
    end
    tweaktau Kv4_fast_GP X
    call Kv4_fast_GP TABFILL X 6000 0
    setfield Kv4_fast_GP X_A->calc_mode {NO_INTERP}
    setfield Kv4_fast_GP X_B->calc_mode {NO_INTERP}

    float tauhmax   = {tauhmax_Kv4f} / {dq10_Kv4}   
    float tauhmin   = {tauhmin_Kv4f} / {dq10_Kv4}
    float Kh        = {Kh_Kv4}
    float V0h       = {V0h_Kv4}
    float Ktauh1    = {Ktauh1_Kv4}
    float Ktauh2    = {Ktauh2_Kv4}
    float hinf, tauh
    call Kv4_fast_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        hinf = 1 / (1 + {exp {({V0h} - x) / {Kh} }})
        tauh = tauhmin + ({tauhmax} - {tauhmin}) / ({exp {({V0h} - x) / {Ktauh1}}} + {exp {({V0h} - x) / {Ktauh2}}}) 
        setfield Kv4_fast_GP Y_A->table[{i}] {tauh}
        setfield Kv4_fast_GP Y_B->table[{i}] {hinf}
        x = x + dx
    end
    tweaktau Kv4_fast_GP Y
    call Kv4_fast_GP TABFILL Y 6000 0
    setfield Kv4_fast_GP Y_A->calc_mode {NO_INTERP}
    setfield Kv4_fast_GP Y_B->calc_mode {NO_INTERP}
end

function make_Kv4_slow_GP
    if (({exists Kv4_slow_GP}))
        return
    end
    create tabchannel Kv4_slow_GP
    float npower = 4
    setfield Kv4_slow_GP Ek {EK} Gbar {G_Kv4_slow_GP} Ik 0 Gk 0\
        Xpower {npower_Kv4} Ypower {hpower_Kv4} Zpower 0
    float Kn        = {Kn_Kv4}
    float V0n       = {V0n_Kv4}
    float taunmax   = {taunmax_Kv4} / {dq10_Kv4}
    float taunmin   = {taunmin_Kv4} / {dq10_Kv4}
    float Ktaun1    = {Ktaun1_Kv4}
    float Ktaun2    = {Ktaun2_Kv4}
    float ninf, taun
    float Vhalfn = {V0n} - ({Kn} * {log {(1 / {pow 0.5 {1/{npower_Kv4}}}) - 1}})
    echo "Kv4s actual Vhalf: " {Vhalfn}
     call Kv4_slow_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        ninf = 1 / (1 + {exp {({V0n} - x) / {Kn} }})
        taun = taunmin + ({taunmax} - {taunmin}) / ({exp {({V0n} - x) / {Ktaun1}}} + {exp {({V0n} - x) / {Ktaun2}}}) 
        setfield Kv4_slow_GP X_A->table[{i}] {taun}
        setfield Kv4_slow_GP X_B->table[{i}] {ninf}
        x = x + dx
    end
    tweaktau Kv4_slow_GP X
    call Kv4_slow_GP TABFILL X 6000 0
    setfield Kv4_slow_GP X_A->calc_mode {NO_INTERP}
    setfield Kv4_slow_GP X_B->calc_mode {NO_INTERP}

    float tauhmax   = {tauhmax_Kv4s} / {dq10_Kv4}   
    float tauhmin   = {tauhmin_Kv4s} / {dq10_Kv4}
    float Kh        = {Kh_Kv4}
    float V0h       = {V0h_Kv4}
    float Ktauh1    = {Ktauh1_Kv4}
    float Ktauh2    = {Ktauh2_Kv4}
    float hinf, tauh
    call Kv4_slow_GP TABCREATE Y {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        hinf = 1 / (1 + {exp {({V0h} - x) / {Kh} }})
        tauh = tauhmin + ({tauhmax} - {tauhmin}) / ({exp {({V0h} - x) / {Ktauh1}}} + {exp {({V0h} - x) / {Ktauh2}}}) 
        setfield Kv4_slow_GP Y_A->table[{i}] {tauh}
        setfield Kv4_slow_GP Y_B->table[{i}] {hinf}
        x = x + dx
    end
    tweaktau Kv4_slow_GP Y
    call Kv4_slow_GP TABFILL Y 6000 0
    setfield Kv4_slow_GP Y_A->calc_mode {NO_INTERP}
    setfield Kv4_slow_GP Y_B->calc_mode {NO_INTERP}
end
