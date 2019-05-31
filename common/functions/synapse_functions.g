// Contains several utility functions for adding synapses to the GP model.
//          J.R. Edgerton, 2008-2009

// ***************************************************************************
// ***************************** HELPER FUNCTIONS ****************************
// ***************************************************************************
//    --> These functions shouldn't need to be called directly by the user.

// ***************************************************************************
// HELPER FUNCTION calc_syn_integ(float gmax, float tau1, float tau2)
//
// Calculate the time integral for a single synaptic event given
//        tau1, tau2, and gmax.
//         J.R. Edgerton, 2006

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
    
// ***************************************************************************
// HELPER FUNCTION setup_ttab(str syntype, str comptname, str tabfile, probR)

/* Create timetable and spikegen objects for a single synapse. Fill the 
    timetable from file {tabfile}, connect messages.

Args:
    1. syntype: type of synapse being connected, such as "STN" or "Striatum"
    2. comptname: name of compartment in which synapse is located
    3. tabfile: name of text file containing timetable data
    4. probR: release probability for timetable (prob_act field)
*/


function setup_ttab(syntype, comptname, tabfile, probR)
    str syntype, comptname, tabfile
    float probR
    str syn_name, this_ttab, this_sgen

    str infostruct = {{syntype} @ "info"}
    if ({probR} == 0)
        echo "WARNING: release probability was set to 0. Changing to 1.0!"
        probR = 1.0
    end

    // determine how many timetables already exist.
    int nsyn = {getfield /{infostruct}/{comptname} num_ttabs}
    // determine how many timetables are supposed to exist.
    int nmax = {getfield /{infostruct}/{comptname} num_syns}
    if (nsyn >= nmax)
        echo "Error. You are trying to add another timetable to the " \
            {syntype} " synapse at " {comptname} ", but the requested number " \
            "of timetables has already been connected."
    else
        nsyn = {nsyn} + 1
        this_ttab = {"ttab" @ {nsyn}}
        this_sgen = {"sgen" @ {nsyn}}
        setfield /{infostruct}/{comptname} num_ttabs {nsyn}

        // set up timetable
        create timetable /{infostruct}/{comptname}/{this_ttab}
        setfield /{infostruct}/{comptname}/{this_ttab} \
            maxtime {rundur} act_val 1.0 method 4 fname {tabfile}
        if ({exists /{infostruct}/{comptname}/{this_ttab} prob_act} == 1)
            setfield /{infostruct}/{comptname}/{this_ttab} prob_act {probR}
        else
            echo "WARNING: timetable object does not have prob_act field."
            echo "Check genesis version."
        end

        call /{infostruct}/{comptname}/{this_ttab} TABFILL    

        // set up spikegen
        create spikegen /{infostruct}/{comptname}/{this_sgen}
        setfield /{infostruct}/{comptname}/{this_sgen} output_amp 1 thresh 0.5

        // pass messages only if timetable was successfully filled.
        if ({getfield /{infostruct}/{comptname}/{this_ttab} allocated} == 1)
            syn_name = {getfield /{infostruct}/{comptname} synname}
            addmsg /{infostruct}/{comptname}/{this_ttab}     \
                /{infostruct}/{comptname}/{this_sgen} INPUT activation
            addmsg /{infostruct}/{comptname}/{this_sgen} \
                   {cellpath}/{comptname}/{syn_name} SPIKE
        else
            // this occurs any time the synaptic rate is 0, or if the tabfile 
            // doesn't exist or is not found. It will NOT crash the simulation.
            echo "Warning: timetable for " {syntype} "  synapse at " \
                 {comptname} " has failed allocation. " \
                "This synapse will not be active."
        end
    end
end


// ***************************************************************************
// HELPER FUNCTION setup_ttab_nofile(str syntype, str comptname, float synrate)
/*
    Create timetable and spikegen objects for a single synapse. Fill the
    timetable internally using the rate parameter synrate. Connect the
    objects appropriately.

Args:
    1. syntype: type of synapse being connected, such as "STN" or "Striatum"
    2. comptname: name of compartment in which synapse is located
    3. synrate: average activity rate (Hz) for this synapse
*/

function setup_ttab_nofile(syntype, comptname, synrate)
    str syntype, comptname
    float synrate
    str syn_name, this_ttab, this_sgen

    str infostruct = {{syntype} @ "info"}
    
    // determine how many timetables already exist.
    int nsyn = {getfield /{infostruct}/{comptname} num_ttabs}
    // determine how many timetables are supposed to exist.
    int nmax = {getfield /{infostruct}/{comptname} num_syns}
    if (nsyn >= nmax)
        echo "Error. You are trying to add another timetable to the " \
            {syntype} " synapse at " {comptname} ", but the requested number " \
            "of timetables has already been connected."
    else
        nsyn = {nsyn} + 1
        this_ttab = {"ttab" @ {nsyn}}
        this_sgen = {"sgen" @ {nsyn}}
        setfield /{infostruct}/{comptname} num_ttabs {nsyn}

        // set up timetable
        create timetable /{infostruct}/{comptname}/{this_ttab}
        if ({synrate} > 0)
            setfield /{infostruct}/{comptname}/{this_ttab}     \
                maxtime {rundur} act_val 1.0 method 2    \
                meth_desc1 {1/{synrate}} meth_desc2 0.0 meth_desc3 3

            call /{infostruct}/{comptname}/{this_ttab} TABFILL    
        end

        // set up spikegen
        create spikegen /{infostruct}/{comptname}/{this_sgen}
        setfield /{infostruct}/{comptname}/{this_sgen} output_amp 1 thresh 0.5

        // pass messages only if timetable was successfully filled.
        if ({getfield /{infostruct}/{comptname}/{this_ttab} allocated} == 1)
            syn_name = {getfield /{infostruct}/{comptname} synname}
            addmsg /{infostruct}/{comptname}/{this_ttab}     \
                /{infostruct}/{comptname}/{this_sgen} INPUT activation
            addmsg /{infostruct}/{comptname}/{this_sgen} \
                   {cellpath}/{comptname}/{syn_name} SPIKE
        else
            // this occurs any time the synaptic rate is 0, or if the tabfile 
            // doesn't exist or is not found. It will NOT crash the simulation.
            echo "Warning: timetable for " {syntype} "  synapse at " \
                 {comptname} " has failed allocation. " \
                "This synapse will not be active."
        end
    end
end


// ***************************************************************************
// HELPER FUNCTION copy_ttab_nofile(str syntype, str comptname, str libobj)
/*
    Copy a template timetable object, create a spikegen object. Connect the
    objects to a synchan.

Args:
    1. syntype: type of synapse being connected, such as "STN" or "Striatum"
    2. comptname: name of compartment in which synapse is located
    3. libobj: name of library timetable object to copy
*/

function copy_ttab_nofile(syntype, comptname, libobj)
    str syntype, comptname, libobj
    str syn_name, this_ttab, this_sgen

    str infostruct = {{syntype} @ "info"}
    
    // determine how many timetables already exist.
    int nsyn = {getfield /{infostruct}/{comptname} num_ttabs}
    // determine how many timetables are supposed to exist.
    int nmax = {getfield /{infostruct}/{comptname} num_syns}
    if (nsyn >= nmax)
        echo "Error. You are trying to add another timetable to the " \
            {syntype} " synapse at " {comptname} ", but the requested number " \
            "of timetables has already been connected."
    else
        nsyn = {nsyn} + 1
        this_ttab = {"ttab" @ {nsyn}}
        this_sgen = {"sgen" @ {nsyn}}
        setfield /{infostruct}/{comptname} num_ttabs {nsyn}
        echo "Compartment " {comptname} " now has " {nsyn} " timetables " \
                "connected to synapse " {syn_name}

        // set up timetable
        copy /library/{libobj} /{infostruct}/{comptname}/{this_ttab}

        // set up spikegen
        create spikegen /{infostruct}/{comptname}/{this_sgen}
        setfield /{infostruct}/{comptname}/{this_sgen} output_amp 1 thresh 0.5

        // pass messages only if timetable was successfully filled.
        if ({getfield /{infostruct}/{comptname}/{this_ttab} allocated} == 1)
            syn_name = {getfield /{infostruct}/{comptname} synname}
            addmsg /{infostruct}/{comptname}/{this_ttab}     \
                /{infostruct}/{comptname}/{this_sgen} INPUT activation
            addmsg /{infostruct}/{comptname}/{this_sgen} \
                   {cellpath}/{comptname}/{syn_name} SPIKE
        else
            // this occurs any time the synaptic rate is 0, or if the tabfile 
            // doesn't exist or is not found. It will NOT crash the simulation.
            echo "Warning: timetable for " {syntype} "  synapse at " \
                 {comptname} " has failed allocation. " \
                "This synapse will not be active."
        end
    end
end



// ***************************************************************************
// ******************************* USER FUNCTIONS ****************************
// ***************************************************************************

//FUNCTION create_syntype_infostruct(str syntype,str comptfile,str wtfile)

/* Read a list of compartments that could potentially receive synapses.
   Create a new neutral element called /{syntype}info, in which there is an 
    entry for each compartment in the list. 
    Each of these elements will be given the following    fields:
        scalefactor: synapse scaling factor
        alloc: 1 = does, 0 = does not have a synapse of type synname.
        synname: name of synaptic conductance (e.g., "AMPA"), initially "none"
 
    At the same time that the compartment names are being read, synaptic weight
    values for those compartments will be read from a separate file ({wtfile}). 
    The weights will be stored in the scalefactor field of the info elements.

    If you just want default scale factors of 1 (e.g. all synapses are the same)
    you can put "none" for the {wtfile} argument.
Input Arguments:
    1. syntype: whatever you want to call this class of inputs
        --> I typically use the presynaptic source ("striatum", "STN", etc.)
    2. comptfile: name of a text file containing compartment names, 1 per line.
    3. wtfile: name of a text file containing synaptic weights, 1 per line.
        --> The ordering of the weight file must match the ordering of the
            compartment file.
*/

function create_syntype_infostruct(syntype, comptfile, wtfile)
    str syntype, comptfile, wtfile
    str comptname
    float synwt

    str infostruct = {{syntype} @ "info"}

    if (! {exists /{infostruct}})
        create neutral /{infostruct}
    end

    // open the ascii text files for reading...
    openfile {comptfile} r
    if ({strcmp {wtfile} "none"} != 0)
        openfile {wtfile} r
        synwt = {readfile {wtfile} -linemode}    // first scale factor from file.
    else
        synwt = 1    // default
    end

    // Get the first compartment name.
    comptname = {readfile {comptfile} -linemode}

    // Loop through the files.
    while (! {eof {comptfile}})

        if ({strcmp {wtfile} "none"} != 0)
            if ({eof {wtfile}})
                echo "Error: not enough scale factors for the number of compts."
                closefile {comptfile}
                closefile {wtfile}
                quit
            end
        end

        // Create a neutral element for this compt, give it the proper fields.
        create neutral /{infostruct}/{comptname}
        addfield /{infostruct}/{comptname} scalefactor
        addfield /{infostruct}/{comptname} num_syns
        addfield /{infostruct}/{comptname} synname
        addfield /{infostruct}/{comptname} num_ttabs
        setfield /{infostruct}/{comptname} scalefactor {synwt} \
            num_syns 0 synname "none" num_ttabs 0

        // Now, get the next compartment name and wt factor.
        comptname = {readfile {comptfile} -linemode}
        if ({strcmp {wtfile} "none"} != 0)
            synwt = {readfile {wtfile} -linemode}
        end
    end
    closefile {comptfile}
    if ({strcmp {wtfile} "none"} != 0)
        if (! {eof {wtfile}})
            echo "Error: number of scaling factors is greater than # compts."
            closefile {wtfile}
            quit
        else
            closefile {wtfile}
        end        
    end
end



// ***************************************************************************
// * FUNCTION add_synchans(str syntype, str fname, str syn_name, str libobj) *

/*
Read a list of compartments from a file. For each compartment, check to see if 
    it already has a synapse named syn_name. If not, copy one from the library.
    To make these synapses active, you must call a separate function to 
    create timetables and connect them.

This function can be called multiple times to add compartment lists from 
    different files. If there is any overlap between the compartment lists, 
    an additional timetable and spikegen will be connected for each entry. 
    So the total number of functional synapses will equal the total number
    of list items read in whether or not they are unique.

Input arguments REQUIRED: 
    1. string with name of synapse type ("striatum", "STN", etc)
    2. string with path & name of file containing the compartment names
    3. string with the name to give these synapses within the compartments
        --> usually you would use "AMPA", but you could separate out your AMPA
            synapses into subgroups this way, like "AMPA_bg" and "AMPA_synch" 
            for background and synchronous synapses, respectively.
    4. string variable with the name of the library synapse object to be
        copied into this compartment. e.g.: "AMPA", "stdpAMPA", etc.

USAGE: 
include <this script>
str syntype = "STN"
str fname = "mySTNcompartments.txt"
str sname = "AMPA_STN"
str libobj = "AMPA"
add_synchans {syntype} {fname} {sname} {libobj}
*/

function add_synchans(syntype, synfname, syn_name, libobj)
    str syntype, synfname, syn_name, libobj

    str infostruct = {{syntype} @ "info"}

    // open file to list compartment names of excitatory synapses
    //    File MUST NOT have any blank lines at the end, or function will fail.
    openfile {synfname} r
    str comptname = {readfile {synfname} -linemode}

    int counter = 0;
    while (! {eof {synfname}})
        if (! {exists /{infostruct}/{comptname}})
            echo "Error: " {comptname} " not found in /" {infostruct}
            closefile {synfname}
            quit
        end
        if ({getfield /{infostruct}/{comptname} num_syns} == 0)
            // Copy a synapse from the library to the compartment    
            copy /library/{libobj} {cellpath}/{comptname}/{syn_name}
            addmsg {cellpath}/{comptname}/{syn_name}    \
                     {cellpath}/{comptname} CHANNEL Gk Ek
            addmsg {cellpath}/{comptname} \
                     {cellpath}/{comptname}/{syn_name} VOLTAGE Vm
            if ({isa stdpSynchan /library/{libobj}})
                addmsg {cellpath}/soma {cellpath}/{comptname}/{syn_name} \
                    VOLTAGE Vm
                setfield {cellpath}/{comptname}/{syn_name} Vspkmsg_idx \
                    {Vspkmsg_stdp}
            end
            // Register the synapse in /{infostruct}
            setfield /{infostruct}/{comptname} synname {syn_name}
            counter = counter + 1
        end
        // Increment num_syns in infostruct
        setfield /{infostruct}/{comptname} num_syns \
            {{getfield /{infostruct}/{comptname} num_syns} + 1}
        // Get next compartment name from file
        comptname = {readfile {synfname} -linemode}
    end
    closefile {synfname}
    return {counter}
end


// ***************************************************************************
// * FUNCTION add_synchans_combofile(str syntype, str fname, str syn_name, str libobj) *

/*
Read a combined compartment names & numbers list (each line has first the
    compartment name, then the # of synapses to go in that compartment).
    Copy the appropriate library synchan object to the compartment, set the
    fields in /{infostruct}/{comptname} to enable timetable setup.
    --> This function does NOT set up or connect timetables. That must be
        handled separately.

Input arguments REQUIRED: 
    1. string with name of synapse type ("striatum", "STN", etc)
    2. string with path & name of file containing the compartment names & #s
    3. string with the name to give these synapses within the compartments
        --> usually you would use "AMPA", but you could separate out your AMPA
            synapses into subgroups this way, like "AMPA_bg" and "AMPA_synch" 
            for background and synchronous synapses, respectively.
    4. string variable with the name of the library synapse object to be
        copied into this compartment. e.g.: "AMPA", "stdpAMPA", etc.

USAGE: 
include <this script>
str syntype = "STN"
str fname = "mySTNcompartments.txt"
str sname = "AMPA_STN"
str libobj = "AMPA"
add_synchans {syntype} {fname} {sname} {libobj}
*/

function add_synchans_combofile(syntype, synfname, syn_name, libobj)
    str syntype, synfname, syn_name, libobj

    str infostruct = {{syntype} @ "info"}

    // open file to list compartment names and numbers
    //    File MUST NOT have any blank lines at the end, or function will fail.
    openfile {synfname} r
    str thisline = {readfile {synfname} -linemode}

    // parse the line into comptname and #syns
    str comptname = {getarg {arglist {thisline}} -arg 1}
    float nsyns = {getarg {arglist {thisline}} -arg 2}

    int counter = 0;
    while (! {eof {synfname}})
        if (! {exists /{infostruct}/{comptname}})
            echo "Error: " {comptname} " not found in /" {infostruct}
            closefile {synfname}
            quit
        end
        if ({getfield /{infostruct}/{comptname} num_syns} != 0)
            echo "Error: " {comptname} " already has synapse of this type.\n"
            closefile {synfname}
            quit
        else
            // Copy a synapse from the library to the compartment    
            copy /library/{libobj} {cellpath}/{comptname}/{syn_name}
            addmsg {cellpath}/{comptname}/{syn_name}    \
                     {cellpath}/{comptname} CHANNEL Gk Ek
            addmsg {cellpath}/{comptname} \
                     {cellpath}/{comptname}/{syn_name} VOLTAGE Vm
            if ({isa stdpSynchan /library/{libobj}})
                addmsg {cellpath}/soma {cellpath}/{comptname}/{syn_name} \
                    VOLTAGE Vm
                setfield {cellpath}/{comptname}/{syn_name} Vspkmsg_idx \
                    {Vspkmsg_stdp}
            end
            // Register the synapse in /{infostruct}
            setfield /{infostruct}/{comptname} synname {syn_name}

            // Set the num_syns field
            setfield /{infostruct}/{comptname} num_syns {nsyns}

            counter = counter + 1
        end
        // Get next compartment name & #syns from file
        thisline = {readfile {synfname} -linemode}
        if (! {eof {synfname}})
            comptname = {getarg {arglist {thisline}} -arg 1}
            nsyns = {getarg {arglist {thisline}} -arg 2}
        end
    end
    closefile {synfname}
    return {counter}
end



// ***************************************************************************
// FUNCTION pad_synchan_list(str syntype, str fname, str syn_name, str libobj, int num_syns)

/*Determine how many synapses of type syntype have already been added, then 
    add additional synapses to compartments listed in {fname} until the total 
    number of synapses = {num_syns}

Input arguments REQUIRED: 
    1. string variable with syntype ("striatum", "STN", etc.)
    2. path & name of file containing the compartment names
    3. string variable with the name to give this group of synchans
        --> e.g. "AMPA_stn", "AMPA_background", etc.
    4. string with name of synchan object to copy from library
    5. integer: total number of synapses desired.
        --> The number added here will be (total - pre-existing)

USAGE: 
include <this script>
str syntype = "STN"
str fname = "myBackgroundSTNlist.txt"
str sname = "AMPA_background"
str libobj = "AMPA"
int num_syns = 100
pad_STN_synlist {syntype} {fname} {sname} {libobj} {num_syns}
*/

function pad_synchan_list(syntype, fname, syn_name, libobj, num_syns)
    str syntype, fname, syn_name, libobj
    int num_syns

    str elname
    str infostruct = {{syntype} @ "info"}

    int syn_count = 0

    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        syn_count = {syn_count} + {getfield {elname} num_syns}
    end

    if ({syn_count} > {num_syns})
        echo "Error: " {num_syns} " synapses were requested, but " {syn_count} \
            " synapses already exist."
        quit
    else
        echo "Number of synapses: " {syn_count}
    end

    echo "Opening file: " {fname}    
    openfile {fname} r
    str comptname = {readfile {fname} -linemode}

    while ({syn_count} < {num_syns})
//        echo "syn_count: " {syn_count} ", num_syns: " {num_syns}
        if ({eof {fname}})
            echo "Error: not enough compartments in list."
            closefile {fname}
            quit
        end
        if ({getfield /{infostruct}/{comptname} num_syns} == 0)
            // Copy a synapse from the library to the compartment    
            echo "Adding " {syn_name} " synapse to compartment: " {comptname}
            copy /library/{libobj} {cellpath}/{comptname}/{syn_name}
            addmsg {cellpath}/{comptname}/{syn_name}    \
                     {cellpath}/{comptname} CHANNEL Gk Ek
            addmsg {cellpath}/{comptname} \
                     {cellpath}/{comptname}/{syn_name} VOLTAGE Vm
            if ({isa stdpSynchan /library/{libobj}})
                addmsg {cellpath}/soma {cellpath}/{comptname}/{syn_name} \
                    VOLTAGE Vm
                setfield {cellpath}/{comptname}/{syn_name} Vspkmsg_idx \
                    {Vspkmsg_stdp}
            end
            // Register the synapse in /{infostruct} and iterate counters
            setfield /{infostruct}/{comptname} synname {syn_name}
        end
        setfield /{infostruct}/{comptname} num_syns \
            {getfield /{infostruct}/{comptname} num_syns} + 1
        syn_count = {syn_count} + 1
        comptname = {readfile {fname} -linemode}
    end
    echo "Total number of " {syntype} " inputs: " {syn_count}
end


// ***************************************************************************
//************* FUNCTION ampscale_syns(str syntype, float Gmax) **************

/* Scale the gmax value of each synapse according to the weight value for that
    synapse's compartment in infostruct. Normalize the scaling such that the 
    mean unitary gmax for the collection of synapses is equal to the original 
    Gmax. This ensures that two different weight distributions will still 
    result in the same total amount of synaptic conductance being applied 
    to the cell, provided that synaptic rates and the number of synapses are 
    the same for both distributions.
    --> Note that if a synapse represents multiple inputs (has multiple
        timetables driving it), all of the inputs will be scaled the same.

Input arguments REQUIRED: 
    1. syntype
    2. Gmax for these synchans

USAGE: 
include <this script>
str syntype = "STN"
float Gmax = 0.25e-9
scale_syns {syntype} {Gmax}
*/

function ampscale_syns(syntype, Gmax)
    str syntype
    float Gmax

    str elname
    int nsyns
    int syn_count = 0
    float totscale = 0
    float meanscale = 0
    str syn_name, comptname

    str infostruct = {{syntype} @ "info"}
    echo "Setting synapse weights..."
    
    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        if ({getfield {elname} num_syns} > 0)
            nsyns = {getfield {elname} num_syns}
            syn_count = {syn_count} + {nsyns}
            totscale = {totscale} + ({getfield {elname} scalefactor}*{nsyns})
        end
    end
    if ({syn_count} == 0)
        echo "No synapses of type " {syntype} " allocated."
        quit
    else
        meanscale = {{totscale} / {syn_count}}
        echo "Number of " {syntype} " syns: " {syn_count}
        echo "Sum total of all scaling factors: "    {totscale} 
        echo "Mean scaling factor: " {meanscale}
        foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
            if ({getfield {elname} num_syns} > 0)
                syn_name = {getfield {elname} synname}
                // get the compartment name without the preceding /{infostruct}/ part.
                comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}
                setfield {cellpath}/{comptname}/{syn_name} gmax        \ 
                    {{Gmax} * {getfield {elname} scalefactor} / {meanscale}}
                echo "synapse " {syn_name} " at " {comptname}    \
                    " has gmax " {getfield {cellpath}/{comptname}/{syn_name} gmax}
            end
        end
    end
end        


// ***************************************************************************
// ***************** FUNCTION scale_synweights(str syntype) *******************

/* Set the synapse.weight value of each synapse equal to the weight value
    in infostruct. Do not normalize the scaling. This is used to initialize
    synapses to the state they were in following a plasticity run such as
    STDP.

Input arguments REQUIRED: 
    1. syntype

USAGE: 
include <this script>
str syntype = "STN"
scale_synweights {syntype}
*/

function scale_synweights(syntype)
    str syntype

    str elname
    int nsyns, i
    str syn_name, comptname

    str infostruct = {{syntype} @ "info"}
    echo "Setting synapse weights..."
    
    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        echo {elname}
        if ({getfield {elname} num_syns} > 0)
            syn_name = {getfield {elname} synname}
            // get the compt name without the preceding /{infostruct}/ part.
            comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}
            nsyns = {getfield {cellpath}/{comptname}/{syn_name} nsynapses}
            if ({nsyns} < {getfield {elname} num_syns})
                echo "Error. scale_synweights must be called AFTER the " \
                    "timetables have been connected, or the synapse weights " \
                    "will not be scaled correctly!"
            end
            for (i=0; i<{nsyns}; i=i+1)
                setfield {cellpath}/{comptname}/{syn_name} \
                    synapse[{i}].weight {getfield {elname} scalefactor}
                echo "synapse " {syn_name} " number " {i} " at " {comptname} \
                 " now has weight " \
                {getfield {cellpath}/{comptname}/{syn_name} synapse[{i}].weight}
            end
        end
    end
end        


// ***************************************************************************
// FUNCTION read_synweights_fromfile(str syn_name, str wtfile, str comptfile, str cellpath)

/* Read synapse weights from an ascii text file, go through the compartments
    listed in a second text file and set each one's weight to the value in the
    file. Do not normalize the scaling. This is used to initialize
    synapses to the state they were in following a plasticity run such as
    STDP.

    *** THIS WORKS ONLY WHERE THE NUMBER OF TIMETABLES PER SYNAPSE = 1! ***

Input arguments REQUIRED: 
    1. syn_name: name of synchan objects to be set, e.g. "AMPA"
    2. wtfile: name of file from which weights will be read.
    3. comptfile: name of file from which compartment names will be read.
    --> obviously the order of the compartment names needs to match the
            order of the synapse weights.
    4. cellpath: path to cell root element (without trailing forward slash)

USAGE: 
include <this script>
str syn_name = "AMPA"
str wtfile = "./mysynweights.asc"
str comptfile = "./mycomptslist.asc"
str cellpath = "/mycell"
read_synweights_fromfile {syn_name} {wtfile} {comptfile} {cellpath}
*/

function read_synweights_fromfile(syn_name, wtfile, comptfile, cellpath)
    str syn_name, wtfile, comptfile, cellpath

    int nsyns, i
    str comptname, elname
    float synwt

    echo "Setting synapse weights..."
    
    // open the ascii text files for reading...
    openfile {comptfile} r
    openfile {wtfile} r

    // Get the first compartment name & synapse weight
    comptname = {readfile {comptfile} -linemode}
    synwt = {readfile {wtfile} -linemode}

    // Loop through the files.
    while (! {eof {comptfile}})

        if ({eof {wtfile}})
            echo "Error: not enough scale factors for the number of compts."
            closefile {comptfile}
            closefile {wtfile}
            quit
        end
        
        nsyns = {getfield {cellpath}/{comptname}/{syn_name} nsynapses}
        for (i=0; i<{nsyns}; i=i+1)
            setfield {cellpath}/{comptname}/{syn_name} \
                synapse[{i}].weight {synwt}
        end
        echo "synapse " {syn_name} " at " {comptname} " has weight " \
             {getfield {cellpath}/{comptname}/{syn_name} synapse[0].weight}

        // Now, get the next compartment name and wt factor.
        comptname = {readfile {comptfile} -linemode}
        synwt = {readfile {wtfile} -linemode}
    end

    closefile {comptfile}
    if (! {eof {wtfile}})
        echo "Error: number of weight factors is greater than # compts."
        closefile {wtfile}
        quit
    else
        closefile {wtfile}
    end        
end
    
// ***************************************************************************
// FUNCTION read_synapse_timetables(str syntype, str comptfile, str ttab_fbase, str ttab_idx_fname, float syn_probR)

/*
Read a list of compartments from {comptfile}. Each of these compartments must
    already have a synapse in it and the synapse must be registered in the 
    information struct for {syntype}.
    A timetable object will be created in /{infostruct}/{comptname}/
    and the table itself will be read in from a file named
    {ttab_fbase} @ {idx}. 

Input Arguments:
    1.    syntype: name for a group of synapses ("striatum", "STN", etc.)
    2.    comptfile: path and name of compartment list.
            e.g. "../../common/GP1/comptlists/gp1_STNinputcomps.asc"
    3.    ttab_fbase: path and filename base for timetable files without the 
            index #s.
            e.g. "../../utilities/timetables/STN/independent/10_HzSTN"
    4.    ttab_idx_fname: path and filename of the timetable indices for this
            compartment list. The first compartment in comptfile would be 
            assigned the timetable with the index number first found in 
            ttab_idx_fname.
            If there are 2 numbers per line in this file, the first is a group
            number and the second a synapse number (for simulations involving
            multiple synchronous groups of inputs).
    5.    syn_probR: release probability to apply to timetables.

    For example, if the first compartment is 'soma', {ttab_fbase} is 
        "../../utilities/timetables/STN/independent/10_HzSTN", and the first 
        row of {ttab_idx_fname} is 25, then the STN synapse in the soma would
        get the timetable read from 
        "../../utilities/timetables/STN/independent/10_HzSTN_25.asc"

This function can be called multiple times to add compartment lists from 
    different files. If a compartment name appears more than once,
    multiple timetables will be connected to the synchan for that compartment.

USAGE: 
include <this script>
str syntype = "STN"
str comptfile = "mySTNcompts.txt"
str ttab_fbase = "../../utilities/timetables/STN/independent/10_HzSTN"
str ttab_idx_fname = "../../common/GP1/ttab_idxs/ttabidx_STN_indep.asc)"
read_synapse_timetables {syntype} {comptfile} {ttab_fbase} {ttab_idx_fname}
*/

function read_synapse_timetables(syntype, comptfile, ttab_fbase, ttab_idx_fname, syn_probR)
    str syntype, comptfile, ttab_fbase, ttab_idx_fname
    float syn_probR

    str comptname, thisline, tabfile, syn_name, basename
    int groupidx, synidx, dset

    echo "Opening compartment list file " {comptfile}
    echo "Opening timetable idx file " {ttab_idx_fname}
        
    openfile {comptfile} r    
    openfile {ttab_idx_fname} r
    
    // get first compartment name
    comptname = {readfile {comptfile} -linemode}

    // get index number(s) for first timetable
    thisline = {readfile {ttab_idx_fname} -linemode}

    // determine how many index numbers were in the line
    int narg = 0
    foreach tstr ({arglist {thisline}})
        // count the number of arguments per line
        narg = narg + 1;
    end

    // parse the line into groupidx and (if applicable) synidx
    if ({narg} == 3)
        dset = {getarg {arglist {thisline}} -arg 1}
        groupidx = {getarg {arglist {thisline}} -arg 2}
        synidx = {getarg {arglist {thisline}} -arg 3}
        basename = "_dset_" @ {dset} @ "_ctxt_" @ {groupidx} @ \ 
                        "_idx_" @ {synidx} @ ".asc"
    elif ({narg} == 2)
        groupidx = {getarg {arglist {thisline}} -arg 1}
        synidx = {getarg {arglist {thisline}} -arg 2}
        basename = "_groupidx_" @ {groupidx} @ "_synidx_" @ {synidx} @ ".asc"

    elif ({narg} == 1)
        groupidx = {getarg {arglist {thisline}} -arg 1}
        synidx = 0
        basename = "_" @ {groupidx} @ ".asc"
    else
        echo "Error: " {ttab_idx_fname} " should have 1 or 2 arguments per line."
    end

    /* Cycle through the list of compartments, create a timetable and spikegen
        for each, fill the tables, pass messages. */
    while (! {eof {comptfile}})
        if ({eof {ttab_idx_fname}})
            echo "Error: not enough timetable index rows for the # of compartments"
            closefile {comptfile}
            closefile {ttab_idx_fname}
            quit
        end

        // Concatenate strings to form complete timetable file name
        tabfile = {ttab_fbase} @ {basename}
        echo "timetable for synapse at " {comptname} " read from: "
        echo "   " {tabfile}

        // call helper function to create timetable & spikegen
        setup_ttab {syntype} {comptname} {tabfile} {syn_probR}

        // get next compartment name and index numbers
        comptname = {readfile {comptfile} -linemode}

        thisline = {readfile {ttab_idx_fname} -linemode}

        if (!{eof {ttab_idx_fname}})
            if ({narg} == 3)
                dset = {getarg {arglist {thisline}} -arg 1}
                groupidx = {getarg {arglist {thisline}} -arg 2}
                synidx = {getarg {arglist {thisline}} -arg 3}
                basename = "_dset_" @ {dset} @ "_ctxt_" @ {groupidx} @ \ 
                                "_idx_" @ {synidx} @ ".asc"
            elif ({narg} == 2)
                groupidx = {getarg {arglist {thisline}} -arg 1}
                synidx = {getarg {arglist {thisline}} -arg 2}
                basename = "_groupidx_" @ {groupidx} @ "_synidx_" @ {synidx} @ ".asc"
            else
                groupidx = {getarg {arglist {thisline}} -arg 1}
                synidx = 0
                basename = "_" @ {groupidx} @ ".asc"
            end
        end
    end
end


// ***************************************************************************
// FUNCTION pad_timetables_fromfiles(str syntype,str ttab_fbase,str ttab_idx_fname, float syn_probR)

/*
Find all synapses that are registered in /{infostruct} but for which there is no
    timetable and spikegen. Create a timetable and spikegen in the same way 
    that read_synapse_timetables would do.

Input Arguments:
    1. syntype
    2.    ttab_fbase: path and filename base for timetable files without the index #s.
            e.g. "../../utilities/timetables/STN/independent/10_HzSTN"
    3.    ttab_idx_fname: path and filename of the timetable indices for this
            compartment list. The first compartment in comptfile would be 
            assigned the timetable with the index number first found in 
            ttab_idx_fname.
            If there are 2 numbers per line in this file, the first is a group
            number and the second a synapse number (for simulations involving
            multiple synchronous groups of inputs).
    4.    startidx: line number in ttabidx file to start at.
    5.    syn_probR: release probability for timetables.
*/

function pad_timetables_fromfiles(syntype, ttab_fbase, ttab_idx_fname, startidx, syn_probR)
    str syntype, ttab_fbase, ttab_idx_fname
    int startidx
    float syn_probR

    if ({startidx} == 0)
        echo "pad_timetables_fromfiles WARNING:"
        echo "Starting line number for ttab_idx file reads is 0."
        echo "Changing to 1."
        startidx = 1
    end

    if ({syn_probR} == 0)
        echo "pad_timetables_fromfiles WARNING:"
        echo "Synapse release probability is set to 0!"
    end

    str comptname, thisline, tabfile, syn_name, basename, elname
    int groupidx, synidx, dset, i, nsyns, nttabs
    
    str infostruct = {{syntype} @ "info"}

    echo "Opening timetable idx file " {ttab_idx_fname}
        
    openfile {ttab_idx_fname} r

    // jump to starting row
    for (i=1; i<={startidx}; i=i+1)    
        // get index number(s) for first timetable
        thisline = {readfile {ttab_idx_fname} -linemode}
    end

    // determine how many index numbers were in the line
    int narg = 0
    foreach tstr ({arglist {thisline}})
        // count the number of arguments per line
        narg = narg + 1;
    end
    echo "Number of args per line in ttab_idx file: " {narg}

    /* Cycle through the compartments, create a timetable and spikegen for 
        each, fill the tables, pass messages. */
    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        // get the compartment name without the preceding /{infostruct}/ part.
        comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}

        // determine if timetables are needed.
        nsyns = {getfield {elname} num_syns}
        nttabs = {getfield {elname} num_ttabs}

        // create & connect the necessary number of new timetables & spkgens
        while ({nttabs} < {nsyns})
            // Make sure {thisline} is valid
            if ({eof {ttab_idx_fname}})
                echo "Error: not enough ttab index rows for the # of compts"
                closefile {ttab_idx_fname}
                quit
            end

            // parse {thisline} into groupidx and (if applicable) synidx
            if ({narg} == 3)
                dset = {getarg {arglist {thisline}} -arg 1}
                groupidx = {getarg {arglist {thisline}} -arg 2}
                synidx = {getarg {arglist {thisline}} -arg 3}
                basename = "_dset_" @ {dset} @ "_ctxt_" @ {groupidx} @ \ 
                                "_idx_" @ {synidx} @ ".asc"

            elif ({narg} == 2)
                groupidx = {getarg {arglist {thisline}} -arg 1}
                synidx = {getarg {arglist {thisline}} -arg 2}
                basename= "_groupidx_" @ {groupidx} @ "_synidx_" @ {synidx} @ ".asc"
            else
                groupidx = {getarg {arglist {thisline}} -arg 1}
                synidx = 0
                basename = "_" @ {groupidx} @ ".asc"
            end

            // Concatenate strings to form complete timetable file name
            tabfile = {ttab_fbase} @ {basename}
            echo "timetable for synapse at " {comptname} " read from: "
            echo "        " {tabfile}

            // call helper function to create timetable & spikegen
            setup_ttab {syntype} {comptname} {tabfile} {syn_probR}

            // increment counter
            nttabs = {nttabs} + 1

            // get index number(s) for next timetable
            thisline = {readfile {ttab_idx_fname} -linemode}    
        end
    end
    closefile {ttab_idx_fname}
end


// ***************************************************************************
// ******** FUNCTION pad_timetables_nofile(str syntype, float synrate) ********

/* 
Scan through all synapses in /{infostruct}, find those that need a 
    timetable and spikegen, create the objects, fill the timetable using the
    rate parameter, and connect the objects to the synapse.

Input Arguments:
    1. syntype: name of synapse to search for (e.g. "striatum" or "STN")
    2. synrate: mean activation rate for each synapse (in Hz)

    /* Cycle through the compartments, create a timetable and spikegen for 
        each provided they don't already exist, fill the tables, pass messages. */
function pad_timetables_nofile(syntype, synrate)
    str syntype
    float synrate

    str infostruct = {{syntype} @ "info"}
    str elname, comptname
    int nsyns, nttabs

    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        // get the compt name without the preceding /{infostruct}/ part.
        comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}

        // determine if timetables are needed.
        nsyns = {getfield {elname} num_syns}
        nttabs = {getfield {elname} num_ttabs}

        // create & connect the necessary number of new timetables & spkgens
        while ({nttabs} < {nsyns})
            // call helper function to create timetable & spikegen
            setup_ttab_nofile {syntype} {comptname} {synrate}
            echo "New timetable added to " {syntype} " synapse at " {comptname}
            nttabs = {nttabs} + 1
        end
    end
end

// ***************************************************************************
// ** FUNCTION pad_timetables_nofile_multisyns(str syntype, float synrate) ***

/* 
Scan through all synapses in /{infostruct}, find those that need a 
    timetable and spikegen, create the objects, fill the timetable using the
    rate parameter, and connect the objects to the synapse.

Input Arguments:
    1. syntype: name of synapse to search for (e.g. "striatum" or "STN")
    2. synrate: mean activation rate for each synapse (in Hz)

    /* Cycle through the compartments, create a timetable and spikegen for 
        each provided they don't already exist, fill the tables, pass messages. */
function pad_timetables_nofile(syntype, synrate)
    str syntype
    float synrate

    str infostruct = {{syntype} @ "info"}
    str elname, comptname
    int nsyns, nttabs

    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        // get the compt name without the preceding /{infostruct}/ part.
        comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}

        // determine if timetables are needed.
        nsyns = {getfield {elname} num_syns}
        nttabs = {getfield {elname} num_ttabs}

        // create & connect the necessary number of new timetables & spkgens
        while ({nttabs} < {nsyns})
            // call helper function to create timetable & spikegen
            setup_ttab_nofile {syntype} {comptname} {synrate}
            echo "New timetable added to " {syntype} " synapse at " {comptname}
            nttabs = {nttabs} + 1
        end
    end
end

// ***************************************************************************
// FUNCTION pad_timetables_nofile_synch(syntype, synrate, idx)

/* 
Scan through all synapses in /{infostruct}, find those that need a 
    timetable and spikegen. Create a SINGLE template timetable and fill its
    table using the synrate parameter. All of the synapses will share this
    common timetable, so they will be synchronous. For each needy synapse,
    copy the template timetable, create a spikegen, connect the objects.

Input Arguments:
    1. syntype: name of synapse to search for (e.g. "striatum" or "STN")
    2. synrate: mean activation rate for each synapse (in Hz)
    3. idx: index number give template ttab, in case there are multiple of
            these.
*/

function pad_timetables_nofile_synch(syntype, synrate, idx)
    str syntype
    float synrate
    int idx, nsyns, nttabs

    str infostruct = {{syntype} @ "info"}
    str elname, comptname

    str libobj = {{syntype} @ "_ttab_" @ {idx}}
    // set up template timetable in library
    create timetable /library/{libobj}
    if ({synrate} > 0)
        setfield /library/{libobj}     \
            maxtime {rundur} act_val 1.0 method 2    \
            meth_desc1 {1/{synrate}} meth_desc2 0.0 meth_desc3 3

        call /library/{libobj} TABFILL    
    end
    
    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        // get the compt name without the preceding /{infostruct}/ part.
        comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}

        // determine if timetables are needed.
        nsyns = {getfield {elname} num_syns}
        nttabs = {getfield {elname} num_ttabs}

        // create & connect the necessary number of new timetables & spkgens
        while ({nttabs} < {nsyns})
            // call helper function to create timetable & spikegen
            copy_ttab_nofile {syntype} {comptname} {libobj}
            echo "Library timetable copied to synapse at " {comptname}
            nttabs = {nttabs} + 1
        end
    end
end


// ***************************************************************************
// FUNCTION pad_timetables_ratescale(str syntype, float synrate, float scl)

/* 
For each compartment receiving an input, determine the ratio of that compt's
    size to the mean size of the whole group of compartments.
    Configure the timetable to use input rate = {synrate} * {size_ratio}
    
Scan through all synapses in /{infostruct}, find those that need a 
    timetable and spikegen, create the objects, fill the timetable using the
    rate parameter, and connect the objects to the synapse.

Input Arguments:
    1. syntype: name of synapse to search for (e.g. "striatum" or "STN")
    2. synrate: mean activation rate for each synapse (in Hz)
    3. scl: rate scaling strategy selector: 
        --if 2, normalize by compartment length
        --if 3, multiply by the num_syns field in /{infostruct}/{comptname}
        --otherwise, normalize by compartment surface area (default)

    /* Cycle through the compartments, create a timetable and spikegen for 
        each provided they don't already exist, fill the tables, pass messages. */
function pad_timetables_ratescale(syntype, synrate, scl)
    str syntype
    float synrate, scl

    str infostruct = {{syntype} @ "info"}
    str elname, comptname
    float d, l, surf
    float totsize, meansize, thisrate, size_ratio
    int num_compts, nsyns, nttabs

    if ({scl} != 3)
        totsize = 0.0
        num_compts = 0
        // cycle through input compartments, determine mean size
        foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
            addfield {elname} surfarea
            addfield {elname} len
            // get the compt name without the preceding /{infostruct}/ part.
            comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}
            //echo {comptname}
            d = {getfield {cellpath}/{comptname} dia}
            if ({strcmp {comptname} "soma"} == 0)
                l = d;
            else
                l = {getfield {cellpath}/{comptname} len}
            end
            surf = {PI}*{l}*{d}
            if ({scl} == 2)
                totsize = {totsize} + {l}
            else
                totsize = {totsize} + {surf}
            end
            num_compts = {num_compts} + 1
            setfield {elname} surfarea {surf}
            setfield {elname} len {l}
        end

        meansize = {totsize} / {num_compts}
        if ({scl} == 2)
            echo "Normalizing by compartment length..."
            echo "Total length (microns): " {{totsize}*1e6}
            echo "Mean length (microns): " {{meansize}*1e6}
        else
            echo "Normalizing by compartment surface area..."
            echo "Total input surface area (microns2): " {{totsize}*1e12}
            echo "Mean S.A. of input compts (microns2): " {{meansize}*1e12}
        end
        echo "Total number of " {syntype} " input events per second: " \
            {{synrate}*{num_compts}}
    end

    foreach elname ({el /{infostruct}/##[OBJECT=neutral]})
        //echo {elname}
        if ({scl} == 2)
            size_ratio = {getfield {elname} len}/{meansize}
        elif ({scl} == 3)
            size_ratio = {getfield {elname} num_syns}
        else
            size_ratio = {getfield {elname} surfarea}/{meansize}
        end
        thisrate = {{synrate}*{size_ratio}}
        // get the compt name without the preceding /{infostruct}/ part.
        comptname = {strsub {elname} {"/" @ {infostruct} @ "/"} ""}
        //echo {comptname}
        echo "Compartment " {comptname} " receives " {thisrate} \
            " " {syntype} " input events per second." 

        if ({scl} == 3)
            nsyns = 1
            nttabs = 0
        else
            // determine how many timetables are needed.
            nsyns = {getfield {elname} num_syns}
            nttabs = {getfield {elname} num_ttabs}
        end

        // create & connect the necessary number of new timetables & spkgens
        while ({nttabs} < {nsyns})
            // call helper function to create timetable & spikegen
            setup_ttab_nofile {syntype} {comptname} {thisrate}
            echo "Rate scaled timetable added to synapse " \
                {syntype} " at " {comptname}
            nttabs = {nttabs} + 1
        end
    end
end


// ***************************************************************************
// FUNCTION add_syns_const(str syntype, str fname, str syn_name, float syn_rate, float Echan, float Gsyn, float tauRise, float tauFall)

/*
Read a list of compartments from a file. For each one, add a leak
    channel object with reversal potential Echan and with constant conductance
    equal to the average conductance that would be expected given the
    rate, unitary size & time-course of the synapse being replaced.

Input arguments REQUIRED: 
    1. string with synapse type ("STN", "striatum", etc.)
    2. string with path & name of file containing the compartment names
    3. string with the name to give these synapses within the compartments
    4. float with the average activation rate (in Hz) for the synapse being
        replaced
    5. float with the reversal potential (in volts) of the synapse.
    6. float with the unitary Gmax of the synapse being replaced.
    7. float with the rise tau (synchan tau1) of the synapse being replaced.
    8. float with the decay tau (synchan tau2) of the synapse being replaced.

USAGE: 
include <this script>
str syntype = "STN"
str fname = "STNinputcompts.asc"
str sname = "AMPA_const"
float synrate = 10
float Echan = 0.0
float Gmax = 0.25e-09
float t1 = 0.001
float t2 = 0.003

add_syns_const {syntype} {fname} {sname}    \
                         {synrate} {Echan} {Gmax} {t1} {t2}
*/

function add_syns_const(syntype, fname, s_name, syn_rate, Echan, Gsyn, t1, t2)
    str syntype, fname, s_name
    float syn_rate, Echan, Gsyn, t1, t2

    str compt
    float gsyn

    str infostruct = {{syntype} @ "info"}

    // calculate average synaptic conductance for given rate
    gsyn = {syn_rate} * { calc_syn_integ {Gsyn} {t1} {t2} }

    echo "Implementing synaptic conductance as tonic leak..."
    echo "Conductance to be added to each input compt (nS): " {gsyn}

    openfile {fname} r
    compt = {readfile {fname} -linemode}

    while (! {eof {fname}})
        if (! {exists /{infostruct}/{compt}})
            echo "Error: " {compt} " not found in /" {infostruct}
            closefile {fname}
            quit
        end
        if ({getfield /{infostruct}/{compt} alloc} == 0)
            // copy leak channel from library to compartment, set values.
            copy /library/leakchan {cellpath}/{compt}/{s_name}
            setfield {cellpath}/{compt}/{s_name} \
                Gbar {gsyn} Ek {Echan}

            // pass messages
            addmsg {cellpath}/{compt} {cellpath}/{compt}/{s_name} VOLTAGE Vm
            addmsg {cellpath}/{compt}/{s_name} {cellpath}/{compt} CHANNEL Gk Ek    
                
            // diagnostic
            echo "Constant synapse added to compt: " {compt}

            // Register the synapse in /{infostruct}
            setfield /{infostruct}/{compt} alloc 1 synname {s_name}

            // Get next compartment name
            compt = {readfile {fname} -linemode}
        end
    end
    closefile {fname}
end


// ***************************************************************************
// FUNCTION add_syns_const_ratescaled(str syntype, str fname, str syn_name, float syn_rate, float Echan, float Gsyn, float tauRise, float tauFall)

/*
Read a list of compartments from a file. Go through the compartments and
    compute the total surface area of all. Then for each one, add a leak
    channel object with reversal potential Echan and with constant conductance
    equal to the average conductance that would be expected given the synapse
    rate, the compartment's surface area relative to that of the others, and
    the unitary size & time-course of the synapse being replaced.

Input arguments REQUIRED: 
    1. string with synapse type ("STN", "striatum", etc.)
    2. string with path & name of file containing the compartment names
    3. string with the name to give these synapses within the compartments
    4. float with the average activation rate (in Hz) for the synapse being
        replaced
    5. float with the reversal potential (in volts) of the synapse.
    6. float with the unitary Gmax of the synapse being replaced.
    7. float with the rise tau (synchan tau1) of the synapse being replaced.
    8. float with the decay tau (synchan tau2) of the synapse being replaced.

USAGE: 
include <this script>
str syntype = "striatum"
str fname = "myGABAcompartments.txt"
str sname = "GABA_const"
float synrate = 10
float Echan = -0.08
float Gmax = 1e-09
float t1 = 0.003
float t2 = 0.012

add_syns_const_ratescaled {syntype} {fname} {sname}    \
                         {synrate} {Echan} {Gmax} {t1} {t2}
*/

function add_syns_const_ratescaled(syntype, fname, s_name, syn_rate, Echan, Gsyn, t1, t2)
    str syntype, fname, s_name
    float syn_rate, Echan, Gsyn, t1, t2

    int i, num_compts
    str compt
    float d,l,surf
    float thisrate, gsyn, gsyn_uni
    float totsurf, meansurf

    str infostruct = {{syntype} @ "info"}

    echo "check1"

    // get sum total surface area of all compartments receiving input.
    totsurf = 0
    openfile {fname} r
    compt = {readfile {fname} -linemode}
    while (! {eof {fname}})
        num_compts = {num_compts} + 1
        d = {getfield {cellpath}/{compt} dia}
        if ({strcmp {compt} "soma"} == 0)
            l = d;
        else
            l = {getfield {cellpath}/{compt} len}
        end
        surf = {PI}*{d}*{l}
        totsurf = {totsurf} + {surf}
        compt = {readfile {fname} -linemode}
    end
    closefile {fname}

    meansurf = {totsurf}/{num_compts}

    // calculate average synaptic conductance for rate of 1 Hz
    gsyn_uni = { calc_syn_integ {Gsyn} {t1} {t2} }

    echo "Implementing synaptic conductance as tonic leak..."
    echo "Total input surface area (microns2): " {{totsurf}*1e12}
    echo "Mean surface area of input compts (microns2): " {{meansurf}*1e12}
    echo "Total number of input events per second: " {{syn_rate}*{num_compts}}
    echo "Events per second for average compartment: " {syn_rate}
    echo "Events per second per square micron: " {{syn_rate}*{num_compts}/({totsurf}*1e12)}
    echo "Total amount of tonic conductance to be applied (nS): " {1e9*{gsyn_uni}*{syn_rate}*{num_compts}}
    echo "Tonic conductance density (S/m2): " {{gsyn_uni}*{syn_rate}*{num_compts}/{totsurf}}

    //cycle through each compartment and add correct amount of leakchan
    openfile {fname} r
    for (i = 1; i <= {num_compts}; i = i + 1)
        compt = {readfile {fname} -linemode}
        if (! {exists /{infostruct}/{compt}})
            echo "Error: " {compt} " not found in /" {infostruct}
            closefile {fname}
            quit
        end

        if ({getfield /{infostruct}/{compt} alloc} == 0)
            d = {getfield {cellpath}/{compt} dia}
            if ({strcmp {compt} "soma"} == 0)
                l = d;
            else
                l = {getfield {cellpath}/{compt} len}
            end
            surf = {PI}*{d}*{l}

            // input rate for this compartment based on surface area
            thisrate = {syn_rate} * {surf} / {meansurf}
         
            // conductance to add is product of unitary average and rate:
            gsyn = {{thisrate}*{gsyn_uni}}
    
            // copy leak channel from library to compartment, set values.
            copy /library/leakchan {cellpath}/{compt}/{s_name}
            setfield {cellpath}/{compt}/{s_name} \
                Gbar {gsyn} Ek {Echan}

            // pass messages
            addmsg {cellpath}/{compt} {cellpath}/{compt}/{s_name} VOLTAGE Vm
            addmsg {cellpath}/{compt}/{s_name} {cellpath}/{compt} CHANNEL Gk Ek    
                
            // Register the synapse in /{infostruct}
            setfield /{infostruct}/{compt} alloc 1 synname {s_name}

            // diagnostic

            echo "COMPT: " {compt} 
            echo "Surface Area (um^2): " {{surf}*1e12}
            echo "Normalized Input Rate (Hz): " {thisrate}
            echo "Constant conductance added (nS): " {{gsyn}*1e9}
        end
    end
    closefile {fname}
end

