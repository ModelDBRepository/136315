// GENESIS SIMULATION IMPLEMENTATION FILE

// *********************** CONSTRUCT LIBRARY & CELL *************************
// Now that all params have been established, create library objects.
//    Intrinsic params should be left alone from this point forward.
include ../common/biophysics/GPsyns.g
include ../common/library/GPcomps.g
include ../common/library/GP_library.g
make_GP_library {mtype}

//create morphology using .p file and compartment objects from library
readcell {morph_fname} {cellpath} -hsolve

// *************************** SIMULATION CLOCKS ****************************
//set up clocks
setclock 0 {dt}          // simulation
setclock 1 {1e-5}     // output--> dictates both output traces and resolution of
                    //    event detectors such as spike_history
setclock 2 {1e-4}    // use for streaming data outputs

/*
// *************************** SYNAPTIC INPUTS ****************************
//add synapses to appropriate compartments    
include ../common/functions/synapse_functions.g

// *** AMPA INPUT ***
    create_syntype_infostruct "STN" {dendfilename} "none"

    add_synchans "STN" {dendfilename} "AMPA_STN" "AMPA"

    // re-seed random number generator
    randseed {rseed_STN}

    // setup timetables & connect to syns
    if ({excscaling} == 0)
        pad_timetables_nofile "STN" {STN_rate}
    elif ({excscaling} == 1)
        // scale rate in proportion to compartment surface area
        pad_timetables_ratescale "STN" {STN_rate} 1
    elif ({excscaling} == 2)
        // scale rate in proportion to compartment length
        pad_timetables_ratescale "STN" {STN_rate} 2
    else
        echo "Error... Unknown value for excscaling parameter."
        quit
    end

// *** GABA INPUT ***
    create_syntype_infostruct "Striatum" {dendfilename} "none"

    add_synchans "Striatum" {dendfilename} "GABA_Str" "GABA"

    // re-seed random number generator
    randseed {rseed_Str}

    // setup timetables & connect to syns
    if ({inhscaling} == 0)
        pad_timetables_nofile "Striatum" {striatum_rate}
    elif ({inhscaling} == 1)
        // scale rate in proportion to compartment surface area
        pad_timetables_ratescale "Striatum" {striatum_rate} 1
    elif ({inhscaling} == 2)
        // scale rate in proportion to compartment length
        pad_timetables_ratescale "Striatum" {striatum_rate} 2
    else
        echo "Error... Unknown value for inhscaling parameter."
        quit
    end
//            ************** END SYNAPSE SETUP **************
*/

// ***************************  CURRENT INJECTION ***************************
create pulsegen /pulse
setfield /pulse                      \
        level1      {{cip_pA}*1e-12} \
        width1      1                \
        delay1      2                \
        delay2      50               \
        baselevel   0                \
        trig_mode   0

addmsg /pulse {cellpath}/soma INJECT output
// ************************* END CURRENT INJECTION **************************


// *****************  IMPLEMENT CHANNEL DENSITY SCALING...  *****************
include ../common/functions/utility_functions.g
if ({chanscale_select} == 2)
    scale_chandens_exp {dendfilename} {chanscale_fname} {cellpath} \
                        "Na_fast_GP" 1.0 {scalemin} {scaletau}
else
    echo "No channel scaling. Conductance densities will stay as default."
    // don't do any scaling--> uniform conductance density.
end

// ****************************** EVENT OUTPUTS *******************************
// NOTE: These outputs must be set up before the hines solver block, or else
//    everything will run as normal but you won't get any event outputs.

if ({save_events} == 1 || {save_STN_ttabs} == 1 || {save_Str_ttabs} == 1)
    // make sim-specific folder to hold event time outputs
    str dirname = {{localdir} @ {basefilename}}
    mkdir {dirname}
end

if ({save_events} == 1)
    echo "saving detected spike events to history files..."
    // Create and connect spikehistory objects...
    str savecompts = "../common/GP1/comptlists/GP1_spkinit.asc"
    save_compt_spiketimes {savecompts} {cellpath} -0.03 0.005 {{dirname} @ "/"}
end

// ****************************** HINES SOLVER *******************************
int cmode, smode
if ({is_hsolved} > 0)
    if ({is_hsolved} == 1)
        if ({save_gSTN} == 1 || {save_gStr} == 1)
            cmode = 4
            smode = 2
        else
            cmode = 3
            smode = 0
        end
    else
        cmode = 0
        smode = 0
    end
    //set up hines solver
        //** Make sure all hsolved elements (channels, synapses, etc) have been
        //  created, inserted and set up before calling SETUP for the solver!
    setfield {cellpath}                         \
        path {cellpath}/##[][TYPE=compartment]  \
            comptmode       1           \
            chanmode        {cmode}     \
            calcmode        0           \
            outclock        1           \
            storemode       {smode}
        call {cellpath} SETUP
        setmethod 11
else
    // Not recommended! Use chanmode 0 if necessary.
end

// **************************** STREAMING OUTPUTS *****************************

if ({save_voltage} == 1)
     // Voltage data...
    create disk_out /out_v
    useclock /out_v 2
    setfield /out_v filename {{localdir} @ {filename_v}} flush 0 append 0 leave_open 1
    
    // save Vm from all compartments in outputcompsfname
    openfile {outputcompsfname} r
    readcompartment = {readfile {outputcompsfname}}
    while (! {eof {outputcompsfname}})
        if ({is_hsolved} == 1)
            hstr = {findsolvefield {cellpath} {cellpath}/{readcompartment} Vm}
            addmsg {cellpath} /out_v SAVE {hstr}
        else
            addmsg {cellpath}/{readcompartment} /out_v SAVE Vm
        end
        readcompartment = {readfile {outputcompsfname}}
    end
    closefile {outputcompsfname}
end

/*
if ({save_gSTN} == 1)
     // STN conductance data...
    create disk_out /out_gSTN
    useclock /out_gSTN 2
    setfield /out_gSTN filename {{localdir} @ {filename_gSTN}} flush 0 append 0 leave_open 1
    
    // save STN syn conductances
    openfile {dendfilename} r
    readcompartment = {readfile {dendfilename}}
    str syn
    while (! {eof {dendfilename}})
        if ({exists {cellpath}/{readcompartment}/AMPA})
            syn = "AMPA"
        elif ({exists {cellpath}/{readcompartment}/AMPA_STN})
            syn = "AMPA_STN"
        elif ({exists {cellpath}/{readcompartment}/AMPAsynch})
            syn = "AMPAsynch"
        else
            echo "Expected AMPA synapse not found in " {readcompartment}
            quit
        end
        if ({is_hsolved} == 1)
            hstr={findsolvefield {cellpath} {cellpath}/{readcompartment}/{syn} Gk}
            addmsg {cellpath} /out_gSTN SAVE {hstr}
        else
            addmsg {cellpath}/{readcompartment}/{syn} /out_gSTN SAVE Gk
        end
        readcompartment = {readfile {dendfilename}}
    end
    closefile {dendfilename}
end

if ({save_gStr} == 1)
    // Striatum conductance data...
    create disk_out /out_gStr
    useclock /out_gStr 2
    setfield /out_gStr filename {{localdir} @ {filename_gStr}} flush 0 append 0 leave_open 1
    
    // save striatum syn conductances
    openfile {dendfilename} r
    readcompartment = {readfile {dendfilename}}
    while (! {eof {dendfilename}})
        hstr = {findsolvefield {cellpath} {cellpath}/{readcompartment}/GABA_Str Gk}
        addmsg {cellpath} /out_gStr SAVE {hstr}
        readcompartment = {readfile {dendfilename}}
    end
    closefile {dendfilename}
end
*/

reset

/*
// **************************** TIMETABLE OUTPUTS *****************************
if ({save_STN_ttabs} == 1)
    echo "Writing timetables to files..."
    // Write synaptic timetables to files
    include ../../common/write_ttabs_separatefiles
    write_timetables {dendfilename} {{dirname} @ "/STN"} "/STNinfo" "ttab"
end

if ({save_Str_ttabs} == 1)
    if ({inhtype} != 0)
        write_timetables {dendfilename} {{dirname} @ "/Striatum"} "/Striatuminfo" "ttab"
    end
end
*/

// **************************** RUN SIMULATION ********************************
// Print data file base name--> used by data file retrieval scripts.
echo "done with output setup."
echo {"DATFILES: " @ {basefilename}}
step {rundur} -time
quit
