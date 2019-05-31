//GENESIS utility script for exporting the times of synaptic excitation in the
//    GP neuron model.
// Author: J. Edgerton, 08/2004
// Modified: J. Edgerton, 09/2008, make function take input parameters
// Modified: J. Edgerton, 08/2010, support synapses with multiple timetables
// Each timetable gets its own output file. Each file is 1 column of ascii txt.

/*
Parameters:
    str comptsfile: ascii file list of names of the compartments from which
        timetables will be saved. 
    str outfilebase: base name of new files to be created.
        --> will be appended with _1.asc, _2.asc, ...
    str ttpath: partial path to timetable elements
        --> will be appended with /{comptname}/{ttname}
        --> All timetables should be located in neutral elements named
            /.../{compartment name}/{ttname}
        --> e.g.: /inputs/STN/soma/timetable
        --> e.g.: /STNinfo/p2b2b1b2b1[79]/ttab
    str ttname: name of timetable elements, e.g. "timetable" or "ttab"
*/
    
function write_timetables(comptsfile, outfilebase, ttpath, ttname)
    str comptsfile, outfilebase, ttpath, ttname

    openfile {comptsfile} r
    str fname
    str thiscompt, thisttab
    int n, k, nttabs

    thiscompt = {readfile {comptsfile} -linemode}
    while (! {eof {comptsfile}})
        nttabs = {getfield {ttpath}/{thiscompt} num_ttabs}
        for (k=1; k<={nttabs}; k=k+1)
            thisttab = {ttpath} @ "/" @ {thiscompt} @ "/" @ {ttname} @ {k}
            fname = {outfilebase} @ "_" @ {thiscompt} @ "_" @ {k} @ ".asc"
            openfile {fname} w
            if (! {exists {thisttab}})
                echo "Error in write_timetables."
                echo "Timetable " {thisttab} " not found"
            end
            ce {thisttab}
            for (n = 0; n <= {getfield ./ maxpos} - 1; n = n+1)
                writefile {fname} {getfield ./ *timetable[{n}]}
            end
            // add one more entry, then start new line.
            closefile {fname}
        end
        thiscompt = {readfile {comptsfile} -linemode}
    end
    closefile {comptsfile}
end
    
    
