
"""
Created on Chien-Chen Lai lab, Institute of Molecular Biology, National Chung Hsing University, Taichung, 40227, Taiwan.

@author: Wun-Yu Lin

This python script is used to integrate MS1 features of Full scan and DDA
Please run the script and input parameter in console 

Pararmeter details:
m/z tolerance: a integer, defining the maximal tolerated m/z deviation in parts per million (ppm) in the same MS1 features of Full and DDA
rt tolerance: a integer, defining the maximal tolerated retention time deviation in senconds in the same MS1 features of Full and DDA
DDA intensity threshold: Minimumn required intensity threshold of DDA MS1 features
DDA number required sample: Minimum required number of sample in DDA MS1 features

"""
#______________________________________________________________________________

#call module
import pandas as pd

#read xlse file
Full_file_path = input("Please input absolute path of Fullscan.xlsx:")
DDA_file_path = input("Please input absolute path of DDA.xlsx:")
mz_tol = input("m/z tolerance, example:20 (ppm), please input number:")
mz_tol = float(mz_tol)
rt_tol = input("rt tolerance, example:30 (secends), please input number:")
rt_tol = int(rt_tol)
DDA_low_intensity = input("DDA intensity threshold, please input number:")
DDA_low_intensity = int(DDA_low_intensity)
DDA_low_intensity_num = input("DDA number required sample, please input number:")
DDA_low_intensity_num = int(DDA_low_intensity_num)
Fullscan_data = pd.read_excel(Full_file_path)
DDA_data = pd.read_excel(DDA_file_path)
print(Full_file_path)
print(DDA_file_path)


DDA_rep3 = pd.DataFrame()
DDA_rep1 = pd.DataFrame()
DDA_data_index = list (DDA_data.index)
##Classify MS1 feature of DDA
for a in DDA_data_index:
    A= list(DDA_data.iloc[a,2:])
    B= len([i for i in A if i>= DDA_low_intensity])
    if B >= DDA_low_intensity_num:
        DDA_rep3 = DDA_rep3.append(DDA_data.iloc[a:a+1])
    else:
        DDA_rep1 = DDA_rep1.append(DDA_data.iloc[a:a+1]) 
print("No of total DDA MS1 features:",len(DDA_data.index))
print("No. of DDA MS1 features after filter",len(DDA_rep3.index))

#Match Fullscan and DDA features
mz1 = Fullscan_data["mz"]
mz2 = DDA_rep3["mz"]
rt1 = Fullscan_data["rt"]
rt2 = DDA_rep3["rt"]
Full_DDA_match = pd.DataFrame()
Full_DDA_match_index =[]
DDA_rep3_index = list (DDA_rep3.index)
Full_index = list (Fullscan_data.index)

print("Starting match Fullscan and DDA features")
for x in DDA_rep3_index:
    mzdiff = abs((mz1-mz2[x])/mz2[x]*1000000)
    rtdiff = abs(rt1-rt2[x])
    for y in Full_index:
        if mzdiff[y] <= mz_tol and rtdiff[y] <= rt_tol:
            Full_DDA_match_index = Full_DDA_match_index+[x]

Full_DDA_match_index = set(Full_DDA_match_index)
Full_DDA_match = Full_DDA_match.append(DDA_rep3.loc[Full_DDA_match_index])
Full_DDA_No_match = DDA_rep3.drop(Full_DDA_match_index)


para_rows = [["m/z tolerance",mz_tol],["rt tolerance",rt_tol],["DDA intensity threshold",DDA_low_intensity],["DDA number required sample",DDA_low_intensity_num],["Fullscan file path",Full_file_path],["DDA file path",DDA_file_path]]
para_df = pd.DataFrame(para_rows, columns = ["parameters","value"])

with pd.ExcelWriter('Fullscan_DDA_integration_results.xlsx') as writer:  
    Full_DDA_match.to_excel(writer, sheet_name='Full_DDA_overlap')
    Full_DDA_No_match.to_excel(writer, sheet_name='Full_DDA_No_overlap')
    DDA_rep3.to_excel(writer, sheet_name='DDA MS1 features after filter')
    para_df.to_excel(writer, sheet_name='parameters')
print("finish")
