set view equal xy
set size sq
set palette
set palette color 
set view map
set style data pm3d
set style function pm3d
#set palette defined (-0.0010"blue",0"white",0.0020"red")
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model CMY
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model XYZ 
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model YIQ
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV 
#set palette rgbformulae 3, 2, 2

#set palette 
#set size sq
#set palette color
#set view map
#set style data pm3d
#set style function pm3d
#set palette defined (0"black",0.5"red")
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model CMY
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model XYZ 
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model YIQ
#set palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV 
#set palette rgbformulae 3, 2, 2
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi)

unset ztics
set ztics format ""
set xtics 10 
set xtics nomirror
set ytics nomirror
set ytics 10
set label 'residue' at -20,40
set label 'residue' at 40,-20 rotate by 90
set tics out  
unset key 
set view 180,90

set colorbox  vertical # LINHA QUE FAZ A LEGENDA
set cbtics 0.0005
set cbrange[0.000:0.001]

set term png

set out 'plot_INTER-CMAP-B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.png'  # FICHEIRO DE OUTPUT

splot [:98][:98]'INTER-CMAP-B2M_D76N-I1-pH7p2---B2M-D76N-I2-pH7p2_2000_T-0.15_n-clashes_0.5_energy_hydropathy+electrostatics+hydrogen-bonds_new.dat' matrix w pm3d   # MUDAR FICHEIRO DE INPUT
