mol new membrane_water.psf
mol addfile quater.dcd waitfor -1

set filename "membrane.txt"

#mol new ionized.pdb type pdb waitfor all
# set sel [atomselect top "((resid >=1 and resid <=28) or (resid >=134 and resid <=140)) and name P"]
set sel [atomselect top "resname DPYF DPYR DPYL and z > 0.0"]
# set indices [$sel get index]

$sel writepsf membrane.psf

set sel [atomselect top "mass > 0.9 and mass < 1.1"]
$sel set name H

set sel [atomselect top "mass > 11.7 and mass < 12.2"]
$sel set name C

set sel [atomselect top "mass > 13.8 and mass < 14.2"]
$sel set name N

set sel [atomselect top "resname DPYF DPYR DPYL and z > 0.0"]
$sel writexyz membrane.xyz

set file [open $filename.ind w]
foreach indices [$sel get index] {
 puts $file "$indices"
}

# flush $file
close $file
