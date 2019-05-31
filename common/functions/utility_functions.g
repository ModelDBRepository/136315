// ***************************************************************************
// FUNCTION scale_chandens_exp(str compt_fname, str scale_fname, str cellpath, str channame, float scalemax, float scalemin, float scaletau)

/*
Read compartment names from file {compt_fname}, read x-axis data from file
    {scale_fname}, generate scaling factors with as exponential decay function
    of x-axis: y = {scalemin} + [{scalemax} * exp(-{x} / {scaletau})]. 
    Multiply channel density by the scaling factor for each compartment.
*/

function scale_chandens_exp(compt_fname, scale_fname, cellpath, channame, scalemax, scalemin, scaletau)
    str compt_fname, scale_fname, cellpath, channame
    float scalemax, scalemin, scaletau

    str thiscompt
    float thisX, thisweight, prevG, newG, GtotPrev, GtotNew

    openfile {compt_fname} r
    openfile {scale_fname} r
    
    thiscompt = {readfile {compt_fname} -linemode}
    thisX = {readfile {scale_fname} -linemode}
    thisweight = {scalemin} + ({scalemax} * {exp {-{thisX}/{scaletau}} })
    GtotNew = 0
    GtotPrev = 0

    while (! {eof {compt_fname}})
        if ({eof {scale_fname}})
            echo "Error in scale_chan_density: not enough scale factors."
            quit
        end
        prevG = {getfield {cellpath}/{thiscompt}/{channame} Gbar}
        echo {cellpath}"/"{thiscompt}"/"{channame}" previous: " {prevG}
        echo {cellpath}"/"{thiscompt}"/"{channame}" scale factor: " {thisweight}
        GtotPrev = {GtotPrev} + {prevG}

        setfield {cellpath}/{thiscompt}/{channame} Gbar {{thisweight} * {prevG}}

        newG = {getfield {cellpath}/{thiscompt}/{channame} Gbar}
        echo {cellpath}"/"{thiscompt}"/"{channame}" new: " {newG}
        GtotNew = {GtotNew} + {newG}        

        thiscompt = {readfile {compt_fname} -linemode}
        thisX = {readfile {scale_fname} -linemode}
        thisweight = {scalemin} + ({scalemax} * {exp {-{thisX}/{scaletau}} })
    end
    echo "Previous Gtot: " {GtotPrev}
    echo "New Gtot: " {GtotNew}
    closefile {compt_fname}
    closefile {scale_fname}
end

// ***************************************************************************
// FUNCTION scale_chandens_lin(str compt_fname, str scale_fname, str cellpath, str channame, float y0, float y_at_xmax)

/*
Read compartment names from file {compt_fname}, read x-axis data from file
    {scale_fname}, generate scaling factors with as linear function
    of x-axis: y = {y0} + [ {slope} * {x} ] 
    Multiply channel density by the scaling factor for each compartment.
    --> Use slope is calculated from y0 & y_at_xmax:
        --> y0 is the scaling factor at x == 0 (x-intercept)
        --> y_at_xmax is y at the largest value of x in the scale_fname file.
        Slope is calculated as (y_at_xmax - y0) / xmax.
        --> calculating slope helps avoid negative scaling factors.
*/

function scale_chandens_lin(compt_fname, scale_fname, cellpath, channame, y0, y_at_xmax)
    str compt_fname, scale_fname, cellpath, channame
    float y0, y_at_xmax

    str thiscompt
    float thisX, thisweight, prevG, newG

    float fileMax = 0

    openfile {scale_fname} r
    thisX = {readfile {scale_fname} -linemode}
    while (! {eof {scale_fname}})
        if ({thisX} > {fileMax})
            fileMax = {thisX}
        end
        thisX = {readfile {scale_fname} -linemode}
    end
    closefile {scale_fname}
    echo "Max Xvalue in file: " {fileMax}
    
    float slope = ({y_at_xmax} - {y0}) / {fileMax}
    echo "Linear Decay Slope: " {slope} " per x-axis unit."

    openfile {compt_fname} r
    openfile {scale_fname} r
    
    thiscompt = {readfile {compt_fname} -linemode}
    thisX = {readfile {scale_fname} -linemode}
    thisweight = {y0} + ({slope} * {thisX})

    while (! {eof {compt_fname}})
        if ({eof {scale_fname}})
            echo "Error in scale_chan_density: not enough scale factors."
            quit
        end
        prevG = {getfield {cellpath}/{thiscompt}/{channame} Gbar}
        echo {cellpath}"/"{thiscompt}"/"{channame}" previous: " {prevG}
        echo {cellpath}"/"{thiscompt}"/"{channame}" scale factor: " {thisweight}

        setfield {cellpath}/{thiscompt}/{channame} Gbar {{thisweight} * {prevG}}

        newG = {getfield {cellpath}/{thiscompt}/{channame} Gbar}
        echo {cellpath}"/"{thiscompt}"/"{channame}" new: " {newG}
        
        thiscompt = {readfile {compt_fname} -linemode}
        thisX = {readfile {scale_fname} -linemode}
        thisweight = {y0} + ({slope} * {thisX})
    end
    closefile {compt_fname}
    closefile {scale_fname}
end


// ***************************************************************************
// FUNCTION save_compt_spiketimes(str compt_fname, str cellpath, float threshold, float refract, str outfilebase)

/*
Read compartment names from file {compt_fname}, create and connect spikegen and
    eventHistory elements to record voltage spike times.

ARGS:
    str compt_fname: name of file listing compartment names
    str cellpath: path to cell root
    float threshold: voltage (in volts) that triggers a spike event
    float refract: minimum time (seconds) after one spike event is registered
        before another can be registered (make sure refract > max spike width
        at threshold).
    str outfilebase: file name to which event times are written will be
        {outfilebase}{thiscompt}_spkhist.asc
*/

function save_compt_spiketimes(compt_fname, cellpath, threshold, refract, outfilebase)
    str compt_fname, cellpath
    float thresh, refract
    str outfilebase

    str thiscompt

    if (! {exists /events})
        create neutral /events
    end    

    openfile {compt_fname} r
    
    thiscompt = {readfile {compt_fname} -linemode}

    while (! {eof {compt_fname}})
        if (! {exists /events/{thiscompt}})
            create neutral /events/{thiscompt}
        end
        create spikegen /events/{thiscompt}/spkdetect
        setfield /events/{thiscompt}/spkdetect output_amp 1 \
             thresh {threshold} abs_refract {refract}
        addmsg {cellpath}/{thiscompt} /events/{thiscompt}/spkdetect INPUT Vm

        create spikehistory /events/{thiscompt}/spkevents
        str fname = {{outfilebase} @ {thiscompt} @ "_spkhist.asc"}
        setfield /events/{thiscompt}/spkevents filename {fname} initialize 1 \
            leave_open 1 flush 1 ident_toggle 1
        addmsg /events/{thiscompt}/spkdetect /events/{thiscompt}/spkevents \
            SPIKESAVE

        thiscompt = {readfile {compt_fname} -linemode}
    end
    closefile {compt_fname}
end


// **************************************************************************
//    FUNCTION randseed_fromfile(str fname, int line_num)
//     Open file fname, get value from line_num, use it to reseed urandom.

function randseed_fromfile(fname, line_num)
    str fname
    int line_num

    openfile {fname} r
    int thisline = {readfile {fname} -linemode}
    int currpos = 1

    while (currpos < line_num)
        if ({eof {fname}})
            echo "Error in randseed_fromfile(): end of file reached."
            currpos = line_num
        else
            thisline = {readfile {fname} -linemode}
            currpos = {currpos} + 1
        end
    end
    
    closefile {fname}
    echo "Reseeding random number generator with " {thisline}
    randseed {thisline}
end
        

// ***************************************************************************
// FUNCTION add_cip_uniform(str compt_fname, str cellpath, float idens, float istart, float idur)

/*
Read compartment names from file {compt_fname}, insert a current pulse generator
    into each compartment within compt_fname with the amplitude scaled by the
    compartment surface area to achieve the specified current density
    --> current density = idens (pA/um2)
    --> start time = istart (s)
    --> duration = idur (s)
*/

function add_cip_uniform(compt_fname, cellpath, idens, istart, idur)
    str compt_fname, cellpath
    float idens, istart, idur

    str thiscompt
    float d, l, surf, icompt, itot

    // Create neutral structure to hold pulsegen objects.
    if ({exists /cip})
        echo "Error in add_cip_uniform. Structure /cip already exists."
        quit
    else
        create neutral /cip
    end

    openfile {compt_fname} r
    
    thiscompt = {readfile {compt_fname} -linemode}

    while (! {eof {compt_fname}})
        d = {getfield {cellpath}/{comptname} dia}
        if ({strcmp {comptname} "soma"} == 0)
            l = d;
        else
            l = {getfield {cellpath}/{comptname} len}
        end
        surf = 1e12 * {PI}*{l}*{d}    // surface area in square microns
        echo "Compartment " {thiscompt} " has surface area " {surf} " um2."
        icompt = {idens}*{surf}        // pA to be injected
        echo "Compartment " {thiscompt} " will receive " {icompt} " pA."
        
        create pulsegen /cip/{thiscompt}
        setfield /cip/{thiscompt}             \
            level1        {1e-12*{icompt}}    \
            width1        {idur}                \
            delay1        {istart}            \
            baselevel    0                    \
            trig_mode    0
        addmsg /cip/{thiscompt} {cellpath}/{thiscompt} INJECT output
        thiscompt = {readfile {compt_fname} -linemode}
    end
    closefile {compt_fname}
end
    
