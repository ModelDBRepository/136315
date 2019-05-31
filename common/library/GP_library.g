//make_GP_library: create the library of components for GP simulation

function make_GP_library(model_select)
    int model_select

    if (!{exists /library})
        create neutral /library
        disable /library
    end

    pushe /library

        make_Na_fast_GP
        make_Na_slow_GP
        make_Kv2_GP
        make_Kv3_GP
        make_Kv4_fast_GP
        make_Kv4_slow_GP
        make_KCNQ_GP
        make_Ca_GP_conc
        make_Ca_GP_nernst
        make_Ca_HVA_GP
        make_SK_GP
        make_h_HCN_GP
        make_h_HCN2_GP
        make_GP_syns
        if ({model_select} == 2)
            echo "making library for NaK GP models (NaF and Kv2 only)"
            echo "Making GP compartments now..."
            make_GP_comps_NaK
        else
            echo "making library for full GP models (all 9 channels)"
            echo "Making GP compartments now..."
            make_GP_comps
        end
    pope
end

