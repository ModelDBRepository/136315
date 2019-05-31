// GENESIS SETUP FILE

silent

//load variable parameter values from environment variable
if ({getenv GENESIS_PAR_ROW} == "")
  echo "*********************************************************************"
  echo "Error: This script needs to read the parameters from the environment "
  echo "        variable GENESIS_PAR_ROW. Set the variable prior to running"
  echo "        the script. Aborting simulation."
  echo "*********************************************************************"
  quit
end

str parrow = {getenv GENESIS_PAR_ROW}

echo "Parameter row: " {parrow}

int i
str tstr, hstr, readcompartment

// ************************  RUN-SPECIFIC PARAMS  **************************
//    Put all parameter settings that will change from run to run here.

// MODEL PROPERTIES
    // mtype: 1 for 9-channel model, 2 for 2-channel (NaK) model
    int mtype = {getarg {arglist {parrow}} -arg 1}

    // INITIALIZE BIOPHYSICS VARIABLES. CAN OVERWRITE WITH PARAMS.
    if ({mtype} == 2)
        include ../common/biophysics/GP1_active_NaK.g
    else
        include ../common/biophysics/GP1_active.g
    end
    include ../common/GP1_constants.g
    include ../common/biophysics/GP1_passive.g
    include ../common/biophysics/GPchannels.g
    include ../common/biophysics/GPsyns_params.g

    /* COMMENT
    ALL intrinsic model params have now been initialized and set. 
    They can be safely overwritten any time between now and the calling of
    the make_GP_library function. Once the library has been created, 
    parameter values are set and cannot be changed except with 
    explicit calls to setfield.
    */

    // chanscale_select: 1 for uniform dendritic gNaF, 2 for declining gradient
    int chanscale_select =    {getarg {arglist {parrow}} -arg 2}

    // scalemin: minimum dendritic gNaF in gradient models is
    //    {initial dendritic gnaf} * {scalemin}
    float scalemin = 0.01    // Gradients would decline to 1% of somatic gNaF if
                            //    the dendrite went far enough.

    // scaletau: gradient length constant (determines how quickly gNaF falls off
    //    with distance from the soma. Units = microns.
    float scaletau =     {getarg {arglist {parrow}} -arg 3}

// CURRENT INJECTION INPUT
    float cip_pA =        {getarg {arglist {parrow}} -arg 4}
    // amplitude of somatic current injection, in pA


// OUTPUT
    // compartments to store Vm from:
    //str outputcompsfname = "../common/comptlists/GP1allcompnames.asc"
    //str outputcompsfname =  "../common/comptlists/soma_only.asc"
    //str outputcompsfname = "../common/comptlists/GP1_spkinit.asc"
    str outputcompsfname = "../common/comptlists/GP1_outputcomps_20060213.asc"

    //filename based on variable parameters
    str basefilename =  {mtype}               @   "_mtype_"       @   \
                        {chanscale_select}    @   "_scaleMeth_"   @   \
                        {scaletau}            @   "_sclTau_"      @   \
                        {cip_pA}              @   "_pAinjected_"  @   \
                        "slice_example_run"

    // Local folder into which data will be written
    str localdir = "data_slice/"

    str filename_v = {basefilename} @ "_v.bin"
//    str filename_gSTN = {basefilename} @ "_gAMPA.bin"
//    str filename_gStr = {basefilename} @ "_gGABA.bin"

    // What will be saved (switches):
    int save_events     = 0    // spike events
    int save_voltage    = 1    // voltage traces
    int save_gSTN       = 0    // All excitatory conductances 
    int save_gStr       = 0    // All inhibitory conductances
    int save_STN_ttabs  = 0    // ASCII text lists of excitatory events
    int save_Str_ttabs  = 0    // ASCII text lists of inhibitory events


// SIMULATION TIMING
    float dt = 1e-5        // time step (in seconds)
    float rundur = 5       // simulation length (in seconds)

// NOW CALL SCRIPT TO HANDLE THE DETAILS
    include setup_slice_example.g
