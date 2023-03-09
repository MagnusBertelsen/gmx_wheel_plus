from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from operator import itemgetter
import seaborn as sns
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.pyplot as mpl
import math
from sklearn.metrics import r2_score



### Parse command line arguments ###
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("-m", "--mode", default="full", help="process = prosess the trajectory for analysis (output = center.xtc), analysis = analyse processed trajectory (center.xtc) to calculate angle distrubution and sidechain pertusion (output = results.txt), wheel = generate helical wheel projection based on results file (default = results.txt, can be changed using -r), analysis+wheel = do analysis + wheel projection, full =  do process + analysis + wheel projection")
parser.add_argument("-f", "--trajectory_name", default="traj.xtc", help="Name of the simulation trajectory")
parser.add_argument("-s", "--tpr_name", default='topol.tpr', help="Name of production simluation tpr-file")
parser.add_argument("-b", "--time_start", default=0, type=int, help="Time for start of analysis in ps / time for once the peptide has reached full insertion into the bilayer in ps")
parser.add_argument("-e", "--end_time", default=500000, type=int, help="Time for end of analysis / end of simultion in ps")
parser.add_argument("-hs", "--helix_start", default=1, type=int, help="First residue used to find the central helical axis. Make sure this resiude is fully part of the helix!")
parser.add_argument("-he", "--helix_end", default=14, type=int, help="Last residue used to find the central helical axis. Make sure this resiude is fully part of the helix!")
parser.add_argument("-r", "--results_input", default='results.txt', help="Input file for wheel projection (only relevant if -m wheel)")
parser.add_argument("-ps", "--projection_start", default="all", help="First resiude shown in the wheel projection")
parser.add_argument("-pe", "--projection_end", default="all", help="Last resiude shown in the wheel projection")


args = vars(parser.parse_args())



### Set up parameters ###

mode=args["mode"]
trajectory_name=args["trajectory_name"]
tpr_name=args["tpr_name"]
time_start=args["time_start"]
end_time=args["end_time"]
helix_start=args["helix_start"]
helix_end=args["helix_end"]
results_input=args["results_input"]
projection_start=args["projection_start"]
projection_end=args["projection_end"]


### Functions ###

def process_trajectory(trajectory_name, tpr_name, time_start, end_time):
  f= open("process_trajectory.sh","w+")
  f.write(
  "#!/usr/bin/bash"+'\n'
  
  ##Arguments
  "trajectory_name=${1}"+'\n'
  "tpr_name=${2}"+'\n'
  "time_start=${3}"+'\n'
  "end_time=${4}"+'\n'
  
  ##Center trajectory on peptide
  "gmx trjconv -f ${1} -s ${2} -o center.xtc -center -pbc mol -ur compact -b ${3} -e ${4} <<EOF"+'\n'
  "Protein"+'\n'
  "System"+'\n'
  "EOF"+'\n')
  
  f.close()  

  os.chmod("process_trajectory.sh", 0o0777)
  
  subprocess.run(["./process_trajectory.sh", str(trajectory_name), str(tpr_name), str(time_start), str(end_time)])
  
  os.remove("process_trajectory.sh")  
  


def get_angle_distrubution(sequence_length, tpr_name, time_start, helix_start, helix_end, end_time):
  f= open("get_angle_distrubution.sh","w+")
  f.write(
  '#!/usr/bin/bash'+'\n'
  
  ##Argument
  'sequence_length=${1}'+'\n'
  'tpr_name=${2}'+'\n'
  'time_start=${3}'+'\n'
  'helix_start=${4}'+'\n'
  'helix_end=${5}'+'\n'
  'end_time=${6}'+'\n'
  
  ##Create index for angle determination
  'echo "del10-30" >> index.txt'+'\n'
  'echo "3 & r${4}-$[${4}+3]" >> index.txt'+'\n'
  'echo "3 & r$[${5}-3]-${5}" >> index.txt'+'\n'
  'for i in $(seq 1 $sequence_length)'+'\n'
  'do'+'\n'
  '	echo "9 & r$i" >> index.txt'+'\n'
  'done'+'\n'
  'echo "del0-9" >> index.txt'+'\n'
  'echo "" >> index.txt'+'\n'
  'echo "q" >> index.txt'+'\n'

  'gmx make_ndx -f ${2} -o wheel.ndx < index.txt'+'\n'
  'rm index.txt'+'\n'

  ##Calculate angle distrubutions for full trajectory
  'for i in $(seq 1 $sequence_length)'+'\n'
  'do'+'\n'
  '	echo "cog of group 0 plus cog of group 1 plus cog of group $[$i+1]" >> gangle_input.txt'+'\n'
  'done'+'\n'

  'gmx gangle -f center.xtc -s ${2} -n wheel.ndx -g1 plane -g2 z -oh anglehis.xvg -oav -xvg none < gangle_input.txt'+'\n'

  ##Calculate angles for rotated helix
  'time=$[${6}-${3}]'+'\n'

  'gmx trjconv -f center.xtc -s ${2} -fit rotxy -dump ${3} -o center_rotxy_start.gro <<EOF'+'\n'
  '3'+'\n'
  '0'+'\n'
  'EOF'+'\n'

  'gmx trjconv -f center.xtc -s ${2} -fit rotxy -dump $[${6}-($time/2)] -o center_rotxy_mid.gro <<EOF'+'\n'
  '3'+'\n'
  '0'+'\n'
  'EOF'+'\n'

  'gmx trjconv -f center.xtc -s ${2} -fit rotxy -dump ${6} -o center_rotxy_end.gro <<EOF'+'\n'
  '3'+'\n'
  '0'+'\n'
  'EOF'+'\n'

  'gmx editconf -f center_rotxy_start.gro -rotate 90 0 0 -o center_rotxy_90_start.gro'+'\n'
  'gmx editconf -f center_rotxy_mid.gro -rotate 90 0 0 -o center_rotxy_90_mid.gro'+'\n'
  'gmx editconf -f center_rotxy_end.gro -rotate 90 0 0 -o center_rotxy_90_end.gro'+'\n'

  'gmx gangle -f center_rotxy_90_start.gro -s center_rotxy_90_start.gro -n wheel.ndx -g1 plane -g2 z -oav angle_rot_start.xvg -xvg none < gangle_input.txt'+'\n'
  'gmx gangle -f center_rotxy_90_mid.gro -s center_rotxy_90_mid.gro -n wheel.ndx -g1 plane -g2 z -oav angle_rot_mid.xvg -xvg none < gangle_input.txt'+'\n'
  'gmx gangle -f center_rotxy_90_end.gro -s center_rotxy_90_end.gro -n wheel.ndx -g1 plane -g2 z -oav angle_rot_end.xvg -xvg none < gangle_input.txt'+'\n'

  'rm gangle_input.txt'+'\n' 
  
  ##Determine N-terminal orientation
  'gmx make_ndx -f ${2} -o terminals.ndx <<EOF'+'\n'
  'del0-40'+'\n'
  'r1'+'\n'
  'r$sequence_length'+'\n'
  'q'+'\n'
  'EOF'+'\n'

  'gmx traj -f center_rotxy_start.gro -s ${2} -n terminals.ndx -x -xvg none -com -ox N-terminal_position.xvg <<<0'+'\n'
  'gmx traj -f center_rotxy_start.gro -s ${2} -n terminals.ndx -x -xvg none -com -ox C-terminal_position.xvg <<<1'+'\n'

  ##Determine insertion leaflet
  'gmx density -f center_rotxy_start.gro -s ${2} -d z -xvg none -o density_system.xvg <<<0'+'\n')


  f.close()  

  os.chmod("get_angle_distrubution.sh", 0o0777)
  
  subprocess.run(["./get_angle_distrubution.sh", str(sequence_length), str(tpr_name), str(time_start), str(helix_start), str(helix_end), str(end_time)])
  
  os.remove("get_angle_distrubution.sh")  



def COG_distance(sequence_length, tpr_name):
  f= open("COG_distance.sh","w+")
  f.write(
  "#!/usr/bin/bash"+'\n'
  
  ##Arguments
  'sequence_length=${1}'+'\n'
  'tpr_name=${2}'+'\n'
  
  'echo "del10-30" >> index_distance.txt'+'\n'
  'for i in $(seq 1 $sequence_length)'+'\n'
  'do'+'\n'
  '	echo "3 & r$i" >> index_distance.txt'+'\n'
  '	echo "9 & r$i" >> index_distance.txt'+'\n'
  'done'+'\n'
  'echo "del1-9" >> index_distance.txt'+'\n'
  'echo "q" >> index_distance.txt'+'\n'

  'gmx make_ndx -f ${2} -o index_distance.ndx < index_distance.txt'+'\n'

  'rm index_distance.txt'+'\n'

  'for i in $(seq 1 2 $[$sequence_length * 2])'+'\n'
  'do'+'\n'
  '	echo "cog of group $i plus cog of group $[$i+1]" >> distance_input.txt'+'\n'
  'done'+'\n'

  'gmx distance -f center.xtc -s ${2} -n index_distance.ndx -oav sidechain_pertusion.xvg -xvg none < distance_input.txt'+'\n'

  'rm distance_input.txt'+'\n')

  f.close()  

  os.chmod("COG_distance.sh", 0o0777)
  
  subprocess.run(["./COG_distance.sh", str(sequence_length), str(tpr_name)])
  
  os.remove("COG_distance.sh")




def pairdist(tpr_name): 
  f= open("pairdist.sh","w+")
  f.write(
  "#!/usr/bin/bash"+'\n'
  
  "tpr_name=${1}"+'\n'
  
  "gmx pairdist -f center.xtc -s ${1} -type max -o pairdist.xvg -xvg none -refgrouping res -selgrouping res <<EOF"+'\n'
  "3"+'\n'
  "9"+'\n'
  "EOF"+'\n')
  
  f.close()  

  os.chmod("pairdist.sh", 0o0777)
  
  subprocess.run(["./pairdist.sh", str(tpr_name)])
  
  os.remove("pairdist.sh")
  
  
  

def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))





def fit_angle(input_file, col):    
    data = pd.read_csv(input_file, delim_whitespace=True, header=None,)
    x_array = np.array(data[0])
    y_array = np.array(data[col])
    
    popt_gauss, pcov_gauss = curve_fit(_1gaussian, x_array, y_array, p0=[0, 90, 40])
    
    fit_sigma = popt_gauss[2]
    fit_cen = popt_gauss[1]
    
    y_pred = _1gaussian(x_array, popt_gauss[0], popt_gauss[1], popt_gauss[2])
    r2 = r2_score(y_array, y_pred)
    
    i=0
    
    while r2 < 0.8:
        popt_gauss, pcov_gauss = curve_fit(_1gaussian, x_array, y_array, p0=[0, i, 20])
        
        fit_sigma = popt_gauss[2]
        fit_cen = popt_gauss[1]
        
        y_pred = _1gaussian(x_array, popt_gauss[0], popt_gauss[1], popt_gauss[2])
        r2 = r2_score(y_array, y_pred)
        
        i = i + 10
        
        if i > 180: 
          fit_cen = 0
          fit_sigma = 0
          r2 = 'NA'
          return (fit_cen, 2*fit_sigma, r2)
          break

    else:
        return (fit_cen, 2*fit_sigma, r2)
    

  
def get_sequence():
  with open('sequence.dat', 'r') as file:
      sequence_raw = file.read().rstrip()

  sequence = []

  for res in sequence_raw:
      sequence.append(res)

  sequence = pd.DataFrame(sequence, columns=['Residue'])
  return sequence
  



def get_sequence_length():
  return len(get_sequence())




def calc_max_distance():
  pairdist = pd.read_csv("pairdist.xvg", delim_whitespace=True, header=None)

  pairdist_selection_list = []
  
  sequence_length = get_sequence_length()

  pairdist_count = 1

  for _ in range(1, sequence_length+1):
     pairdist_selection_list.append(pairdist_count)
     pairdist_count = pairdist_count + sequence_length + 1
   
  pairdist = pairdist[pairdist_selection_list]
  pairdist_mean = pairdist.mean()
  pairdist_mean = pd.DataFrame(pairdist_mean, columns = ['Sidechain length'])
  pairdist_mean = pairdist_mean.reset_index(drop=True)

  return pairdist_mean



def calc_COG_distance():
  distance = pd.read_csv("sidechain_pertusion.xvg", delim_whitespace=True, header=None)
  
  distance_mean = distance.mean()
  distance_mean = distance_mean.drop(axis=0, index = 0)
  distance_mean = pd.DataFrame(distance_mean, columns = ['Sidechain length'])
  distance_mean = distance_mean.reset_index(drop=True)
  
  return distance_mean



def calc_angle():
  angle_dis = []
  sequence_length = get_sequence_length()

  for res in range(1, sequence_length+1):
      angle_dis.append(fit_angle("anglehis.xvg", res))

  angle_dis = pd.DataFrame(angle_dis, columns = ['Angle centre', 'Angle spread (2xSigma)', 'Fit R^2'])

  return angle_dis
  
  
  
def calc_angle_rot():
  angle_rot_start = pd.read_csv("angle_rot_start.xvg", delim_whitespace=True, header=None)
  angle_rot_start = angle_rot_start.transpose()	
  angle_rot_start = angle_rot_start.drop(axis=0, index = 0,)
  angle_rot_start.columns = ['Rotated angle start']  
  angle_rot_start = angle_rot_start.reset_index(drop=True)  
  
  angle_rot_mid = pd.read_csv("angle_rot_mid.xvg", delim_whitespace=True, header=None)
  angle_rot_mid = angle_rot_mid.transpose()	
  angle_rot_mid = angle_rot_mid.drop(axis=0, index = 0,)
  angle_rot_mid.columns = ['Rotated angle mid']  
  angle_rot_mid = angle_rot_mid.reset_index(drop=True)  
  
  angle_rot_end = pd.read_csv("angle_rot_end.xvg", delim_whitespace=True, header=None)
  angle_rot_end = angle_rot_end.transpose()	
  angle_rot_end = angle_rot_end.drop(axis=0, index = 0,)
  angle_rot_end.columns = ['Rotated angle end']  
  angle_rot_end = angle_rot_end.reset_index(drop=True)  
  
  angle_rot = pd.concat([angle_rot_start, angle_rot_mid, angle_rot_end], axis=1)
  angle_rot = angle_rot.mean(axis=1)
  angle_rot = pd.DataFrame(angle_rot, columns = ['Rotated angle'])
  
  return angle_rot
  
  
  
def get_360_angle():
  
  ## X-orientation ##
  N_terminal_position = pd.read_csv("N-terminal_position.xvg", delim_whitespace=True, header=None)
  N_terminal_position = N_terminal_position.iat[0, 1]
  C_terminal_position = pd.read_csv("C-terminal_position.xvg", delim_whitespace=True, header=None)
  C_terminal_position = C_terminal_position.iat[0, 1]
  
  ### Leaflet insertion ###
  density = pd.read_csv('density_system.xvg', delim_whitespace=True, header=None)
  density_length = density[0].max(axis=0)
  index_density_min = np.argmin(density[1])
  density_min = density.iloc[index_density_min] 
  density_min = density_min[0]
  
  angle_dis = calc_angle()
  
  angle_rot = calc_angle_rot()
  
  angle_full = []



  if N_terminal_position > C_terminal_position and density_min < (density_length/2): #works
    for ind in angle_dis.index:
      if angle_rot['Rotated angle'][ind] < 90:
        angle_full.append(360 - angle_dis['Angle centre'][ind] - 180)
      else:
        angle_full.append(angle_dis['Angle centre'][ind] +180)
        #print('A')
      
   
  elif N_terminal_position < C_terminal_position and density_min < (density_length/2): #works
    for ind in angle_dis.index: 
      if angle_rot['Rotated angle'][ind] > 90:
         angle_full.append(360 - angle_dis['Angle centre'][ind] - 180)
      else:
         angle_full.append(angle_dis['Angle centre'][ind] + 180)  
         #print('B')
       
  elif N_terminal_position > C_terminal_position and density_min > (density_length/2): #works
    for ind in angle_dis.index: 
      if angle_rot['Rotated angle'][ind] < 90:
         angle_full.append(360 - angle_dis['Angle centre'][ind])
      else:
         angle_full.append(angle_dis['Angle centre'][ind])
         #print('C')
       
  else: 
    for ind in angle_dis.index: 
      if angle_rot['Rotated angle'][ind] > 90:
       angle_full.append(360 - angle_dis['Angle centre'][ind])
      else:
         angle_full.append(angle_dis['Angle centre'][ind])
         #print('D')   
        
  angle_full = pd.DataFrame(angle_full, columns = ['Angle 360'])
  
  print(angle_full)
  return angle_full
  


def calc_results():
  sequence = get_sequence() 
  angle_dis = calc_angle() 
  angle_rot = calc_angle_rot() 
  angle_full = get_360_angle() 
  distance_max = calc_max_distance()
  
  results = pd.concat([sequence, angle_dis, angle_rot, angle_full, distance_max], axis=1)
  
  results.to_csv('results.txt', index=False, sep=' ')
  
  return results


def get_data(results_input, projection_start, projection_end):
  
  if projection_end == 'all' or projection_start == 'all': 
    data = pd.read_csv(results_input, sep=" ")
    return data
  
  else:
    seqeunce_length = get_sequence_length()
    residue_drop_list =[]
    
    for i in range(1,int(projection_start)):
      residue_drop_list.append(i)
    for i in range(int(projection_end)+1, seqeunce_length+1):
       residue_drop_list.append(i)
    
    data = pd.read_csv(results_input, sep=" ", skiprows=residue_drop_list, header=0)
    return data


def gen_wheel(results_input, projection_start, projection_end):
  
  data = get_data(results_input, projection_start, projection_end)

  #Make graph basic 
  ax = plt.subplot(111, polar=True)
  angles = np.deg2rad(data["Angle 360"])
  heights = data["Sidechain length"]
  widths = np.deg2rad(data["Angle spread (2xSigma)"])
  residues = data["Residue"]

  #Colouring
  label_color_dict = {'H': '#00FFFF',
                      'I': '#b4b4b4',
                      'L': '#b4b4b4',
                      'K': '#00FFFF',
                      'M': '#b4b4b4',
                      'F': '#b4b4b4',
                      'T': '#00FF00',
                      'W': '#b4b4b4',
                      'V': '#b4b4b4',
                      'R': '#00FFFF',
                      'C': '#FFFF00',
                      'Q': '#00FF00',
                      'G': '#FFFF00',
                      'P': '#FFFF00',
                      'Y': '#b4b4b4',
                      'A': '#b4b4b4',
                      'D': '#FF0000',
                      'N': '#00FF00',
                      'E': '#FF0000',
                      'S': '#00FF00',}
                    
  colors = itemgetter(*residues)(label_color_dict)
	
  ## Draw bars
  bars = ax.bar(
      x=angles, 
      height=(heights), 
      width=widths, 
      bottom=0.23,
      linewidth=2, 
      edgecolor="white",
      color=colors,
      alpha=0.5,
      )

  plt.axis('off')

  ## Add labels
  residue_labels = []
  if projection_start == 'all' or projection_end == 'all':
    x = 0
    for i in residues:
      x = x + 1
      residue_labels.append("{}{}".format(x, i))
      
  else:
    x = int(projection_start) -1
    for i in residues:
      x = x + 1
      residue_labels.append("{}{}".format(x, i))    

  labelPadding = 0.15
  for bar, angle, height, label in zip(bars,angles, heights, residue_labels):

      ax.text(
          x=angle, 
          y=bar.get_height() + labelPadding, 
          s=label,  
          va='center',
          ha='center') 

  plt.show()


def main(mode, trajectory_name, tpr_name, time_start, end_time, helix_start, helix_end, results_input, projection_start, projection_end):
 
  if mode == "process":
    process_trajectory(trajectory_name, tpr_name, time_start, end_time)
  
  elif mode == "analysis":
    sequence_length = get_sequence_length()
    get_angle_distrubution(sequence_length, tpr_name, time_start, helix_start, helix_end, end_time)
    pairdist(tpr_name)
    calc_results()
  
  elif mode == "wheel":
    gen_wheel(results_input, projection_start, projection_end)
  
  elif mode == "analysis+wheel":
    sequence_length = get_sequence_length()
    get_angle_distrubution(sequence_length, tpr_name, time_start, helix_start, helix_end, end_time)
    pairdist(tpr_name)
    calc_results()
    gen_wheel(results_input, projection_start, projection_end)
    
  else:
    process_trajectory(trajectory_name, tpr_name, time_start, end_time)
    sequence_length = get_sequence_length()
    get_angle_distrubution(sequence_length, tpr_name, time_start, helix_start, helix_end, end_time)
    pairdist(tpr_name)
    calc_results()
    gen_wheel(results_input, projection_start, projection_end)






if __name__ == '__main__':
  main(mode, trajectory_name, tpr_name, time_start, end_time, helix_start, helix_end, results_input, projection_start, projection_end)   
    


