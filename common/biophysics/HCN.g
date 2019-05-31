//==================================================================
//      H current  (Anomalous rectifier--mixed Na and K current)
//        HCN1/HCN2 heteromeric channels and HCN2 homomeric channels
//        Channel model from Chan et al (2004), J Neurosci 24: 9921-32.
//        Original model from Siegelbaum lab. Wang et al (2002), Neuron 36:
//            451-62. Chen et al (2001), JGP 117: 491-504.
//==================================================================

// --> kinetics for room temperature

// HCN1/2 heteromeric channels:
float mpower_HCN1    = 1
float V0m_HCN1       = -0.0764
float Km_HCN1        = -0.0033
float taumin_HCN1    = 0
float taumax_HCN1    = 14.5    // actual taumax 3.625 with Q10 adjustment
float Ktau1_HCN1     = 0.00656
float Ktau2_HCN1     = -0.00748
float dq10_HCN1      = 4

function make_h_HCN_GP
    if ({exists h_HCN_GP})
               return
    end
    create tabchannel h_HCN_GP
    setfield h_HCN_GP Ek {Eh} Gbar {{G_h_HCN_GP}}  \
        Xpower {mpower_HCN1} Ypower 0 Zpower 0
    float Km         = {Km_HCN1}
    float V0m        = {V0m_HCN1}
    float taumin     = {taumin_HCN1} / {dq10_HCN1}
    float taumax     = {taumax_HCN1} / {dq10_HCN1}
    float Ktau1      = {Ktau1_HCN1}
    float Ktau2      = {Ktau2_HCN1}
    float dq10       = {dq10_HCN1}
    float minf, taum
    call h_HCN_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        minf = 1 / (1 + {exp {({V0m} - x) / {Km} }})
        taum = {taumin} + (({taumax} - {taumin}) / ({exp {({V0m} - x)/{Ktau1}}} + {exp {({V0m}-x)/{Ktau2}}}))
        setfield h_HCN_GP X_A->table[{i}] {taum}
        setfield h_HCN_GP X_B->table[{i}] {minf}
        x = x + dx
    end
    tweaktau h_HCN_GP X
    call h_HCN_GP TABFILL X 6000 0
    setfield h_HCN_GP X_A->calc_mode {NO_INTERP}
    setfield h_HCN_GP X_B->calc_mode {NO_INTERP}
end

// HCN2 homomeric channels:
float mpower_HCN2    = 1
float V0m_HCN2       = -0.0875
float Km_HCN2        = -0.004
float taumin_HCN2    = 0
float taumax_HCN2    = 25.2    // actual taumax 6.3 with Q10 adjustment
float Ktau1_HCN2     = 0.0089
float Ktau2_HCN2     = -0.0082
float dq10_HCN2      = 4

function make_h_HCN2_GP
    if ({exists h_HCN2_GP})
           return
    end
    create tabchannel h_HCN2_GP
    setfield h_HCN2_GP Ek {Eh} Gbar {{G_h_HCN2_GP}}  \
        Xpower {mpower_HCN2} Ypower 0 Zpower 0
    float Km         = {Km_HCN2}
    float V0m        = {V0m_HCN2}
    float taumin     = {taumin_HCN2} / {dq10_HCN2}
    float taumax     = {taumax_HCN2} / {dq10_HCN2}
    float Ktau1      = {Ktau1_HCN2}
    float Ktau2      = {Ktau2_HCN2}
    float dq10       = {dq10_HCN2}
    float minf, taum
    call h_HCN2_GP TABCREATE X {xdivs} {xmin} {xmax}
    x = xmin
    for (i = 0; i <= {xdivs}; i = i + 1)
        minf = 1 / (1 + {exp {({V0m} - x) / {Km} }})
        taum = {taumin} + (({taumax} - {taumin}) / ({exp {({V0m} - x)/{Ktau1}}} + {exp {({V0m} - x)/{Ktau2}}}))
        setfield h_HCN2_GP X_A->table[{i}] {taum}
        setfield h_HCN2_GP X_B->table[{i}] {minf}
        x = x + dx
    end
    tweaktau h_HCN2_GP X
    call h_HCN2_GP TABFILL X 6000 0
    setfield h_HCN2_GP X_A->calc_mode {NO_INTERP}
    setfield h_HCN2_GP X_B->calc_mode {NO_INTERP}
end
