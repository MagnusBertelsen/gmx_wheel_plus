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



**Example output:**
![v3](https://user-images.githubusercontent.com/127429845/224114917-bd2ec12e-78a7-4f62-8d28-ad74ed80cc9e.png)
