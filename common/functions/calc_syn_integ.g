// Script to calculate the time integral for a single synaptic event given
//    tau1, tau2, and gmax.
// J.R. Edgerton, 2006

function calc_syn_integ(gmax, tau1, tau2)
    float gmax, tau1, tau2
    float pktime = ({log {tau1}}-{log {tau2}})*{tau1}*{tau2}/({tau1}-{tau2})
    float pkval = ({gmax} / ({tau1} - {tau2}))    \
        * ({exp {-{pktime} / {tau1}}} - {exp {-{pktime} / {tau2}}})
    float A = gmax/{abs {pkval}}
    // Now have all components of synapse alpha function:
    //    gsyn = ((A * gmax)/(tau1-tau2)) * (exp(-t/tau1) - exp(-t/tau2))

    // Can now integrate alpha function over time. Stop integration when
    //    time reaches 5*tau2, where amplitude < 1% of peak.

    float syninteg = ({A} * {gmax}) / ({tau1} - {tau2})    \
        * (({tau2} * {exp -5}) - ({tau1} * {exp {-5*{tau2}/{tau1}}}) \
        - {tau2} + {tau1})

    return {syninteg}
end
    
