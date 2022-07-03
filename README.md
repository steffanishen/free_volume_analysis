to run:
1. Create the .xyz and index file for dcd
   vmd -dispdev text -eofexit <create_index.tcl> out
2. Carve out the dcd
   catdcd -o membrane.dcd -i membrane.txt.ind -dcd quater.dcd
3. Revise the input.dat
   Note one entry is the number of processors
4. Run:
   sbatch submit.sh
