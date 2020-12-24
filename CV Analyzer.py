# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:13:59 2020

@author: Hannah Teeters
"""
"""This code processes all .txt files in a directory. Before running this code, 
make sure directory contains no outside .txt files. To run the code more than once,
remove "CV Analysis.txt" before beginning. Otherwise an error message will occur""" 

"""Python Tools"""
import os.path
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
df=pd.DataFrame()

"""Input Mass and Molecular Weight here"""
Mass=2
MW=142.00

"""Do not change conversion factor (CF) """
CF= (3.6/96485.3329) 

"""Code must be in same folder as .mpt files"""
get_directory=os.getcwd()
mfile_location=os.path.join(get_directory,'*.mpt')
mfilenames=glob.glob(mfile_location)

"""Empty lists to use throughout loops"""
headings=[]
ac=[]
cc=[]
out_list=[]

"""Opens .mpt files for Python to read"""
for f in mfilenames:    
    file_name=os.path.basename(f)
    split_filename=file_name.split('.')
    file_num=len(split_filename)    
    moutfile=open(f,'r')
    data=moutfile.readlines()
    moutfile.close()
    
    """Formats sweep rate, voltage, and cycle data for use in data frame"""           
    for line in data:
        if 'dE/dt               ' in line:
            sweep_line=line
            data1=sweep_line.split("            ")
            headings.append(data1[1])
        elif 'E1 (V)              ' in line:
            volt1_line=line
            data2=volt1_line.split("             ")
            headings.append(data2[1])
        elif 'E2 (V)              ' in line:
            volt2_line=line
            data3=volt2_line.split("              ")
            headings.append(data3[1])
        elif 'nc cycles           ' in line:
            cycle_line=line
            data4=cycle_line.split("      ")
            headings.append(data4[1])

"""This function moves all of the actual data to a new .txt file free of headers 
but keeps the original files intact"""
def move_lines(original_file, new_file):
    is_skipped = False
    counter=0
    new_file=original_file+'.txt'
    with open(original_file, 'r') as read_obj, open(new_file, 'w') as write_obj:
        for line in read_obj:
            if counter > 53:
                write_obj.write(line)
            else:
                is_skipped = True
            counter=counter+1
    while is_skipped == False:
        os.write(new_file, line)

"""Utilizes function to copy actual data into new .txt files. 
   Headers no longer exist in new files""" 
for f in mfilenames:    
    move_lines(file_name, 'headers.txt')

"""Calls ALL .txt files in directory. """
get_directory=os.getcwd()
tfile_location=os.path.join(get_directory,'*.txt')
tfilenames=glob.glob(tfile_location)

"""Loop to run calculations for ALL .txt files"""
for f in tfilenames:  
    
    """Creates data frame to read files"""
    df=pd.read_table(f, index_col=None, header=None, names=["mode", "ox/red", "error", 
                                             "ctrlchgs", "counter", "Time", 
                                             "ctrlV", "EweV", "Current", 
                                             "Cycle", "Q-Qo/C", "P/W"], delim_whitespace = True)
    
    """Removes unnecessary columns from dataframe"""
    df.drop(["mode", "ox/red", "error", "ctrlchgs", 
             "counter", "ctrlV", "Q-Qo/C", "P/W"], axis=1, inplace=True)
    
    """Converts Cycle and Current values to actual numbers"""
    df.Cycle = pd.to_numeric(df.Cycle, errors='coerce').fillna(0).astype(np.int64)
    df.Current = pd.to_numeric(df.Current, errors='coerce').fillna(0).astype(np.float64)
    
    """Removes Cycles 1&3
    If you want to view cycles 1&3, comment out these lines"""
    df.drop(df[df['Cycle'] < 2].index, inplace = True)
    df.drop(df[df['Cycle'] > 2].index, inplace = True)
    
    """Puts current and time into a separate array"""
    it=df[['Current', 'Time']].to_numpy()   
    
    """Classifies data as either positive or negative based on current"""
    bool_neg= it[:,0] < 0 
    bool_pos= it[:,0] >= 0 
    it_pos=it[bool_pos]
    it_neg=it[bool_neg]
    
    """Calculates Anodic Specific Capacity"""
    yp=it_pos[:,0]
    xp=it_pos[:,1]
    Anod_cap= abs(np.trapz(yp, x=xp/Mass, dx=1, axis=-1))
    
    """Calculates Cathodic Specific Capacity"""
    yn=it_neg[:,0]
    xn=it_neg[:,1]
    Cath_cap= abs(np.trapz(yn, x=xn/Mass, dx=1, axis=-1))
    
    """Calculates remaining values"""    
    CE= (Anod_cap/Cath_cap)
    EsA= (Anod_cap*MW*CF)
    EsC= (Cath_cap*MW*CF)
    ac.append(Anod_cap)
    cc.append(Cath_cap)    
    rca = (Anod_cap/max(ac))
    rcc = (Cath_cap/max(cc))

    """Prints all calculations and plots current vs time
    ~'ro' in plot command shows individual points
    ~round function rounds values to the nearest 10,000th (4th decimal place)"""
    print("Anodic Specific Capacity: ", round(Anod_cap, 4))
    print("Coulombic Efficiency: ", round(CE,4))
    print("E_stored_Anod: ", round(EsA,4))
    print("Rate Capability Anode: ", round(rca,4))
    plt.plot(xp, yp, 'ro')
    plt.show()
    print("Cathodic Specific Capacity: ", round(Cath_cap,4))
    print("E_stored_Cath: ", round(EsC,4))
    print("Rate Capability Cathode: ", round(rcc,4))
    plt.plot(xn, yn, 'ro')
    plt.show()
    
    """Creates data frame and adds values to empty list for every loop"""
    output = pd.DataFrame({'Anodic Spec Cap':[round(Anod_cap, 4)],
                           'Cathodic Spec Cap':[round(Cath_cap,4)],
                           'Coulumbic Efficiency':[round(CE,4)],
                           'E Stored Anodic':[round(EsA,4)],
                           'E Stored Cathodic':[round(EsC,4)],
                           'Rate Cap Anodic': [round(rca,4)],
                           'Rate Cap Cathodic':[round(rcc,4)]})
    out_list.append(output)

""""Transfers list contents to new data frame and numpy array"""
out_df=pd.concat(out_list)
out_dat=out_df.to_numpy()

"""Puts numpy array back into data frame with headers"""
data_frame=pd.DataFrame({'Anodic Specific Capacity (mAh/g)': out_dat[:,0],
                         'Cathodic Specific Capacity (mAh/g)': out_dat[:,1],
                         'Coulumbic Efficiency': out_dat[:,2],
                         '#e- Stored Anode': out_dat[:,3],
                         '#e- Stored Cathode':out_dat[:,4],
                         'Rate Capability Anodic': out_dat[:,5],
                         'Rate Capability Cathodic':out_dat[:,6]})

"""Transfers list of header info into data frame"""
out_headings=pd.DataFrame(headings)

"""Converts data frame to numpy array to reshape"""
out_mat= out_headings.to_numpy()
new_mat=out_mat.reshape((-1, 4))

"""Converts new matrix into data frame"""
header_frame=pd.DataFrame({'Sweep Rate (mV/s)': new_mat[:,0],
                      'Voltage 1 (V)': new_mat[:,1],
                      'Voltage 2 (V)': new_mat[:,2],
                      'Cycles': new_mat[:,3]})
"""Joins header info and calculations in a new data frame"""
new_concat=pd.concat([header_frame, data_frame], axis=1)

"""Writes new data frame to new text file"""
tfile=open('CV Analysis.txt', 'w')
tfile.write(new_concat.to_string())
tfile.close()



    
    
    
        

