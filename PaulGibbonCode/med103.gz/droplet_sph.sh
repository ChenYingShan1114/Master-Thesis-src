#
#  Run directory - change name to suit run description
#
RUNDIR=sphere1
#
if [ -d $RUNDIR ] 
then
  echo "Run directory" $RUNDIR "already exists"
else
  echo "Creating run directory" $RUNDIR
  mkdir $RUNDIR
fi

#
echo 'Running MEDUSA ..'
cd $RUNDIR
# Clean
rm -f *.000 fort.*
touch fort.11 fort.12 fort.13
#
#  Assumes that an input deck 'droplet.ind' exists in the base directory
#
cp ../droplet_sph.ind med.ind
#
#  Execute code
#
../med103 < med.ind 
#
#  Use this form to suppress 'lineprinter' output!
#../med103 < med.ind > med.out
#
#  To monitor run progress type:
#    tail -f $RUNDIR/med.out | grep timestep
#
#  Rename some output files
mv fort.11  med.xrl
#  ion pops
mv fort.12 med.ion
#  hydro snapshots
mv fort.13 med_graph.asc
#  laser intensity and timestep history
mv fort.90 laser.dat
#
echo '...finished'
echo
echo 'Postprocessing ...'
../medgraf
echo '... done'
ls ne* > snap_list  ## Make new snapshot list
echo "Snapshots:"
cat snap_list | while read SNAP 
do
        A=${SNAP##*e}    ## strip 'ne' from string
#        B=`expr $A * 1000`
       echo $A "ps" 
done




