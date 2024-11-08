#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Linux version 3.7
#    	patchlevel 0.1
#    	last modified Thu Mar  4 17:08:06 GMT 1999
#    
#    	Copyright(C) 1986 - 1993, 1998, 1999
#    	Thomas Williams, Colin Kelley and many others
#    
#    	Type `help` to access the on-line reference manual
#    	The gnuplot FAQ is available from
#    		<http://www.uni-karlsruhe.de/~ig25/gnuplot-faq/>
#    
#    	Send comments and requests for help to <info-gnuplot@dartmouth.edu>
#    	Send bugs, suggestions and mods to <bug-gnuplot@dartmouth.edu>
#    
set terminal postscript color 
set output 'ne_snap.ps'
set noclip points
set clip one
set noclip two
set border 31 lt -1 lw 1.000
set xdata
set ydata
set zdata
set x2data
set y2data
set boxwidth
set dummy x,y
set format x "%g"
set format y "%g"
set format x2 "%g"
set format y2 "%g"
set format z "%g"
set angles radians
set nogrid
set key title ""
set key right top Right noreverse box linetype -2 linewidth 1.500 samplen .6 spacing 1 width 0
set nolabel
set noarrow
set nolinestyle
set nologscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
set nopolar
set noparametric
set view 50, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
set nocontour
set clabel '%8.3g'
set mapping cartesian
set hidden3d offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set cntrparam order 4
set cntrparam linear
set cntrparam points 5
set size ratio 0 0.5,0.5
set origin 0,0
set data style lines
set function style lines
set xzeroaxis lt -2 lw 1.000
set x2zeroaxis lt -2 lw 1.000
set yzeroaxis lt -2 lw 1.000
set y2zeroaxis lt -2 lw 1.000
set tics in
set ticslevel 0.5
set ticscale 1 0.5
set mxtics default
set mytics default
set mx2tics default
set my2tics default
set xtics border mirror norotate autofreq 
set ytics border mirror norotate autofreq 
set ztics border nomirror norotate autofreq 
set nox2tics
set noy2tics
set title "" 0.000000,0.000000  ""
set timestamp "" bottom norotate 0.000000,0.000000  ""
set rrange [ * : * ] noreverse nowriteback  # (currently [-0:10] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set xlabel "" 0.000000,0.000000  ""
set x2label "" 0.000000,0.000000  ""
set timefmt "%d/%m/%y\n%H:%M"
set xrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set x2range [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set ylabel "" 0.000000,0.000000  ""
set y2label "" 0.000000,0.000000  ""
set yrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set y2range [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set zlabel "" 0.000000,0.000000  ""
set zrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set zero 1e-08
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale "C"
set size 2,1
set hidden
set multiplot
set size .6,.6
set origin 0.,.5
splot 'rhoe01.snap'
set origin 0.,0.
splot 'rhoe02.snap'
set origin .5,.5
splot 'rhoe03.snap'
set origin .5,0.
splot 'rhoe04.snap'
set nomultiplot

#    EOF
