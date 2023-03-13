# gmx_wheel_plus
gmx_wheel+ is a tool for analysing and visualizing the orientation of helical antimicrobial peptides based on molecular dynamics simulations trajectories generated using the Gromacs simulation suite. 


**Purpose of the tool:**

The general purpose of the tool is to generate more accurate helical wheel "like" projections of antimicrobial peptide based on more advanced _in silico_ data rather than just the primary sequence. 

  

The amphiphilic structure of the helical region of antimicrobial peptides is key to their activity, which has led helical wheel projection to be one of the most simple but effective tools for researches to use in their design process. 
The simple and clear projection of the structure of the peptide allows one to easily analyse the residue composition and select residues for modification when designing new analogs. 

 
With the lowering of the barrier of entry to more advanced techniques such as molecular dynamic simulations, much more advanced _in silico_ data can be generated for an antimicrobial peptide before synthesis and _in vitro_ evaluation. 
However, the data still needs to be analysed and visualized in ways to best facilitate the design process. 




**How it works:**

The gmx_wheel+ tool can be used to analyse the orientation of helical peptides after they have been inserted into a phospholipid bilayer in their surface parallel configuration. The angles between the normal of the bilayer (Z-axis) and the planes created based on the centre of geometry of start of the helix, end of the helix and each amino acid side chain are estimated for each frame using the gmx gangle function. The angle distributions are evaluated using gaussian fits, calculating the mean value and the spread and transformed from 180° to 360° by correcting against three sideways measurements of the angles. In addition, the average maximum distance between the alpha-carbons and the sidechains are calculated for each residue using the gmx pairdist function. The results are graphed using a polar axis bar plot using mathplotlib, with the position of each bar being based on the mean angle for each residue, the width of each bar being based on 2x the spread and the hight is based on the average max protrusion of the sidechain.   



**How to use:**

Copy the gmx_wheel+.py script into a folder containing the following files: 

1) The trajectory file of the peptide-bilayer simulation.  

2) The tpr file for the peptide-bilayer simulation. 

3) A plain text file named sequence.dat, containing only the peptide sequence in single letter code. 

The script can be run through the terminal using different arguments to customise the inputs, functions and outputs. 

Information related to each argument can be seen using the –h flag. 

The scripts have several modes that can be selected using –m, either running the full script or individual parts such as the trajectory processing, analysis or graphing only. 

Names of the input files are specified using the –s and –f flags. Start and end times of the period of the trajectory that should be analysed are specified using –b and –e. Make sure that the peptide has reached an equilibrium position within the bilayer leaflet before starting the analysis. The start and end point of the helical part of the peptide can be specified using the –hs and –he flags. Make sure that these points and the 3 residues after (-hs) and before (-he) are fully helical such that the centre of geometry lines up with the helical axis. Non-ideal helical residues outside of the range are still analysed and shown in the figure if they are within the projection range set by –ps and –pe.  

 

An example of the command line input could look like this: 

>python3 gmx_wheel+_v1.0.py -m full -s md_1_500.tpr -f md_1_500.xtc -b 200000 -e 500000 -hs 2 -he 12 -ps 1 –pe 14  



**Example output:**

Below is an example of the output for the 14-residue antimicrobial peptide Smp14a, inserted into a 1:1 DOPC:DOPG bilayer. 



![v3](https://user-images.githubusercontent.com/127429845/224114917-bd2ec12e-78a7-4f62-8d28-ad74ed80cc9e.png)

Blue resiudes are cationic, green residues are polar and grey resiudes are hydrophobic.
