import os,re,subprocess
import numpy as np
import pandas as pd
import math
import csv


def get_outstreams(log): # gets the compressed stream information at the end of a Gaussian job
    streams = []
    starts,ends = [],[]
    error = "failed or incomplete job" # default unless "normal termination" is in file

    try:
        with open(log+".log") as f:
            loglines = f.readlines()
    except:
        with open(log+".LOG") as f:
            loglines = f.readlines()
    for i in range(len(loglines)):
        if "1\\1\\" in loglines[i]:
            starts.append(i)
        if "@" in loglines[i]:
            ends.append(i)
        if "Normal termination" in loglines[i]:
            error = ""
    if len(starts) != len(ends) or len(starts) == 0: #probably redundant
        error = "failed or incomplete job"
        return(streams,error)
    for i in range(len(starts)):
        tmp = ""
        for j in range(starts[i],ends[i]+1,1):
            tmp = tmp + loglines[j][1:-1]
        #print(tmp)
        streams.append(tmp.split("\\"))
    return(streams,error)


def get_SPE_Pd(filename_nolog, atoms_P, atoms_Pd_Cl2, atoms_R, atoms_backbone):
    SPE_dict = {}
    filename = filename_nolog + '.log'
    atom_keys = ['P1', 'P2', 'Pd', 'Cl1', 'Cl2', 'R1', 'R2', 'R3', 'R4', 'RBack1', 'RBack2']
    atom_values= atoms_P + atoms_Pd_Cl2 + atoms_R + atoms_backbone
    atom_label_dict =dict(zip(atom_keys, atom_values))

# Added 3_10_22 JJD
    HF_E = []
    SPE_dict['SPE_HF_E'] = 'Error'
    stream = get_outstreams(filename_nolog)
    HF_E = []
    for job in stream[0]:
        for line in job:
            pattern = 'HF=-'
            if re.search(pattern, line):
                HF_E.append(line[3:])
    #print(filename_nolog, HF_E)
    SPE_dict['SPE_HF_E'] = HF_E[0]

# Added 3_29_2022 JJD PR1_NBO, R2_NBO, P_Rback_NBO,R3_NBO, R4_NBO, X_Rback_NBO C    4   -0.37018
    R1_NBO = ''
    SPE_dict['R1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R1']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['R1']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    R1_NBO = str.split(line)[2]
                    SPE_dict['R1_NBO'] = R1_NBO
                break

    R2_NBO = ''
    SPE_dict['R2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R2']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['R2']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    R2_NBO = str.split(line)[2]
                    SPE_dict['R2_NBO'] = R2_NBO
                break

    RBack1_NBO = ''
    SPE_dict['RBack1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['RBack1']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['RBack1']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    RBack1_NBO = str.split(line)[2]
                    SPE_dict['RBack1_NBO'] = RBack1_NBO
                break

    R3_NBO = ''
    SPE_dict['R3_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R3']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['R3']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    R3_NBO = str.split(line)[2]
                    SPE_dict['R3_NBO'] = R3_NBO
                break

    R4_NBO = ''
    SPE_dict['R4_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R4']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['R4']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    R4_NBO = str.split(line)[2]
                    SPE_dict['R4_NBO'] = R4_NBO
                break

    RBack2_NBO = ''
    SPE_dict['RBack2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['RBack2']+"\s+-0." +"|"+ "[CON]\s+" + atom_label_dict['RBack2']+"\s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    RBack2_NBO = str.split(line)[2]
                    SPE_dict['RBack2_NBO'] = RBack2_NBO
                break

##31 P NMR
    P1_NMR = ''
    SPE_dict['P1_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P1']+'\s+P\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    P1_NMR = (str.split(line)[4])
                    ##print('P1 NMR: ', P1_NMR)
                    SPE_dict['P1_NMR'] = P1_NMR
                break


    P2_NMR = ''
    SPE_dict['P2_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P2']+'\s+[PN]\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    P2_NMR = (str.split(line)[4])
                    #print('X NMR: ', P2_NMR)
                    SPE_dict['P2_NMR'] = P2_NMR
                break


##Aniosotropic NMR shifts (added by LvD Nov 2021)
    aniso_P1_NMR = ''
    SPE_dict['aniso_P1_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P1']+'\s+P\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    aniso_P1_NMR = (str.split(line)[7])
                    ##print('P1 NMR: ', P1_NMR)
                    SPE_dict['aniso_P1_NMR'] = aniso_P1_NMR
                break


    aniso_P2_NMR = ''
    SPE_dict['aniso_P2_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P2']+'\s+[PN]\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    aniso_P2_NMR = (str.split(line)[7])
                    #print('X NMR: ', P2_NMR)
                    SPE_dict['aniso_P2_NMR'] = aniso_P2_NMR
                break



    HOMO = []
    SPE_dict['Homo']= 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = r"Alpha\s+occ"
            if re.search(pattern, line):
                HOMO = str.split(line)[-1]
                break
        ##print('HOMO:', HOMO)
        SPE_dict['Homo']= HOMO

    LUMO = []
    SPE_dict['Lumo'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = r"Alpha\s+virt\."
            if re.search(pattern, line):
                LUMO = str.split(line)[4]
                break
        ##print('LUMO:', LUMO)
        SPE_dict['Lumo'] = LUMO

    Pd_dipole = []
    SPE_dict['dipole'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = " X= "
            if re.search(pattern, line):
                if str.split(line)[-2] == "Tot=":
                    Pd_dipole = str.split(line)[-1]
                break
        ##print("Pd dipole: ", Pd_dipole)
        SPE_dict['dipole'] = Pd_dipole

    P1_NBO = []
    SPE_dict['P1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            #pattern = "P    " + atom_label_dict['P1']+"    [1-]"
            pattern = "P\s+" + atom_label_dict['P1']+"    [1-]"
            if re.search(pattern, line):
                if str.split(line)[0] == "P":
                    P1_NBO = str.split(line)[2]
                    SPE_dict['P1_NBO'] = P1_NBO
                break
        ##print("P1 NBO:", P1_NBO)


    P2_NBO = []
    SPE_dict['P2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[PN]\s+" + atom_label_dict['P2'] +"    [1-]"
            if re.search(pattern, line):
            ##    if str.split(line)[0] == "[PN]":
                P2_NBO = str.split(line)[2]
                SPE_dict['P2_NBO'] = P2_NBO
                break
        #print("X NBO:", P2_NBO)




    Pd_NBO = []
    SPE_dict['Pd_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            #pattern = "Pd   " + atom_label_dict['Pd']+"   [-0]"
            pattern = "Pd\s+" + atom_label_dict['Pd']+"\s+[-0]"
            if re.search(pattern, line):
                if str.split(line)[0] == "Pd":
                    Pd_NBO = str.split(line)[2]
                    SPE_dict['Pd_NBO'] = Pd_NBO
                break
        ##print("Pd NBO:", Pd_NBO)


    Cl_1_NBO = []                          ############
    SPE_dict['Cl_1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            #pattern = "Cl    " + atom_label_dict['Cl1']+"   [-0]"
            pattern = "Cl\s+" + atom_label_dict['Cl1']+"   [-0]"
            if re.search(pattern, line):
                if str.split(line)[0] == "Cl":
                    Cl_1_NBO = str.split(line)[2]
                    SPE_dict['Cl_1_NBO'] = Cl_1_NBO
                break
        ##print("Cl #1 NBO:", Cl_1_NBO)


    Cl_2_NBO = []
    SPE_dict['Cl_2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            #pattern = "Cl    " + atom_label_dict['Cl2']+"   [-0]"
            pattern = "Cl\s+" + atom_label_dict['Cl2']+"   [-0]"
            if re.search(pattern, line):
                if str.split(line)[0] == "Cl":
                    Cl_2_NBO = str.split(line)[2]
                    SPE_dict['Cl_2_NBO'] = Cl_2_NBO
                break
        ##print("Cl #2 NBO:", Cl_2_NBO)



    P1_Pd_bond_occ = []
    P1_Pd_bond_eng = []
    SPE_dict['P1_Pd_bond_occ'] = 'Error'
    SPE_dict['P1_Pd_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            #pattern = "BD \(\s+1\) P\s+" + atom_label_dict['P1']+"\s+-Pd\s+" + atom_label_dict['Pd']
            pattern = "BD \(\s+1\) P\s+" + atom_label_dict['P1']+"\s+-Pd\s+" + atom_label_dict['Pd']+"|"+"BD \(\s+1\)Pd\s+" + atom_label_dict['Pd'] + " - P\s+" + atom_label_dict['P1']
            if re.search(pattern, line):
                P1_Pd_bond_occ = str.split(line)[8]
                P1_Pd_bond_eng = str.split(line)[9]
                SPE_dict['P1_Pd_bond_occ'] = P1_Pd_bond_occ
                SPE_dict['P1_Pd_bond_eng'] = P1_Pd_bond_eng
                break
        ##print("P1-Pd bond occ :", P1_Pd_bond_occ)
        ##print("P1-Pd bond eng :", P1_Pd_bond_eng)


    P2_Pd_bond_occ = []                      ##########
    P2_Pd_bond_eng = []
    SPE_dict['P2_Pd_bond_occ'] = 'Error'
    SPE_dict['P2_Pd_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) [PN]\s+" + atom_label_dict['P2']+" -Pd\s+" + atom_label_dict['Pd'] +"|"+"BD \(\s+1\)Pd\s+" + atom_label_dict['Pd'] + " - [PN]\s+" + atom_label_dict['P2']
            if re.search(pattern, line):
                P2_Pd_bond_occ = str.split(line)[8]
                P2_Pd_bond_eng = str.split(line)[9]
                SPE_dict['P2_Pd_bond_occ'] = P2_Pd_bond_occ
                SPE_dict['P2_Pd_bond_eng'] = P2_Pd_bond_eng
                break
        #print("X-Pd bond occ :", P2_Pd_bond_occ)
        #print("X-Pd bond eng :", P2_Pd_bond_eng)


    P1_Pd_antibond_occ = [] ### Why are there so many hits for the search phrase, "BD*(   1) P   3 -Pd  74"?
    P1_Pd_antibond_eng = []
    SPE_dict['P1_Pd_antibond_occ'] = 'Error'
    SPE_dict['P1_Pd_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) P\s+" + atom_label_dict['P1'] + " -Pd\s+" + atom_label_dict['Pd'] + "|"+"BD\*\(\s+1\)Pd\s+" + atom_label_dict['Pd'] + " - P\s+" + atom_label_dict['P1']
            if re.search(pattern, line):
                P1_Pd_antibond_occ = str.split(line)[7]
                P1_Pd_antibond_eng = str.split(line)[8]
                SPE_dict['P1_Pd_antibond_occ'] = P1_Pd_antibond_occ
                SPE_dict['P1_Pd_antibond_eng'] = P1_Pd_antibond_eng
                break
        ##print("P1-Pd antibond occ :", P1_Pd_antibond_occ)
        ##print("P1-Pd antibond eng :", P1_Pd_antibond_eng)


    P2_Pd_antibond_occ = []
    P2_Pd_antibond_eng = []
    SPE_dict['P2_Pd_antibond_occ'] = 'Error'
    SPE_dict['P2_Pd_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) [PN]\s+" + atom_label_dict['P2'] + " -Pd\s+" + atom_label_dict['Pd'] + "|"+"BD\*\(\s+1\)Pd\s+" + atom_label_dict['Pd'] + " - [PN]\s+" + atom_label_dict['P2']
            if re.search(pattern, line):
                P2_Pd_antibond_occ = str.split(line)[7]
                P2_Pd_antibond_eng = str.split(line)[8]
                SPE_dict['P2_Pd_antibond_occ'] = P2_Pd_antibond_occ
                SPE_dict['P2_Pd_antibond_eng'] = P2_Pd_antibond_eng
                break
        ##print("X-Pd antibond occ :", P2_Pd_antibond_occ)
        ##print("X-Pd antibond eng :", P2_Pd_antibond_eng)



    Pd_LP_1_occ = [] ### What is the difference between LP 1 and LP 2?
    Pd_LP_1_eng = []### Why don't you include this energy value in your list of parameters?
    SPE_dict['Pd_LP_1_occ'] = 'Error'
    SPE_dict['Pd_LP_1_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP \(\s+1\)Pd\s+" + atom_label_dict['Pd']
            if re.search(pattern, line):
                Pd_LP_1_occ = str.split(line)[5]
                Pd_LP_1_eng = str.split(line)[6]
                SPE_dict['Pd_LP_1_occ'] = Pd_LP_1_occ
                SPE_dict['Pd_LP_1_eng'] = Pd_LP_1_eng
                break
        ##print("Pd LP 1 occ :", Pd_LP_1_occ)
        ##print("Pd LP 1 eng :", Pd_LP_1_eng)


    Pd_LP_2_occ = []
    Pd_LP_2_eng = []
    SPE_dict['Pd_LP_2_occ'] = 'Error'
    SPE_dict['Pd_LP_2_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP \(\s+2\)Pd\s+" + atom_label_dict['Pd']
            if re.search(pattern, line):
                Pd_LP_2_occ = str.split(line)[5]
                Pd_LP_2_eng = str.split(line)[6]
                SPE_dict['Pd_LP_2_occ'] = Pd_LP_2_occ
                SPE_dict['Pd_LP_2_eng'] = Pd_LP_2_eng
                break
        ##print("Pd LP 2 occ :", Pd_LP_2_occ)
        ##print("Pd LP 2 eng :", Pd_LP_2_eng)


    Pd_LP_3_occ = []
    Pd_LP_3_eng = []
    SPE_dict['Pd_LP_3_occ'] = 'Error'
    SPE_dict['Pd_LP_3_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP \(\s+3\)Pd\s+" + atom_label_dict['Pd']
            if re.search(pattern, line):
                Pd_LP_3_occ = str.split(line)[5]
                Pd_LP_3_eng = str.split(line)[6]
                SPE_dict['Pd_LP_3_occ'] = Pd_LP_3_occ
                SPE_dict['Pd_LP_3_eng'] = Pd_LP_3_eng
                break
        ##print("Pd LP 3 occ :", Pd_LP_3_occ)
        ##print("Pd LP 3 eng :", Pd_LP_3_eng)


    Pd_LP_4_occ = []
    Pd_LP_4_eng = []
    SPE_dict['Pd_LP_4_occ'] = 'Error'
    SPE_dict['Pd_LP_4_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP \(\s+4\)Pd\s+" + atom_label_dict['Pd']
            if re.search(pattern, line):
                Pd_LP_4_occ = str.split(line)[5]
                Pd_LP_4_eng = str.split(line)[6]
                SPE_dict['Pd_LP_4_occ'] = Pd_LP_4_occ
                SPE_dict['Pd_LP_4_eng'] = Pd_LP_4_eng
                break
        ##print("Pd LP 4 occ :", Pd_LP_4_occ)
        ##print("Pd LP 4 eng :", Pd_LP_4_eng)


    P1_R1_bond_occ = []
    P1_R1_bond_eng = []
    SPE_dict['P1_R1_bond_occ'] = 'Error'
    SPE_dict['P1_R1_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R1']  + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R1']+" - P\s+"+ atom_label_dict['P1'] + '\s+'## Remember to add spaces or periods to patterns that end with a number
            if re.search(pattern, line):
                P1_R1_bond_occ = str.split(line)[9]
                P1_R1_bond_eng = str.split(line)[10]
                SPE_dict['P1_R1_bond_occ'] = P1_R1_bond_occ
                SPE_dict['P1_R1_bond_eng'] = P1_R1_bond_eng
                break
        ##print("P1-R1 bond occ :", P1_R1_bond_occ)
        ##print("P1-R1 bond eng :", P1_R1_bond_eng)

    P2_R3_bond_occ = []
    P2_R3_bond_eng = []
    SPE_dict['P2_R3_bond_occ'] = 'Error'
    SPE_dict['P2_R3_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R3'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R3']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                P2_R3_bond_occ = str.split(line)[9]
                P2_R3_bond_eng = str.split(line)[10]
                SPE_dict['P2_R3_bond_occ'] = P2_R3_bond_occ
                SPE_dict['P2_R3_bond_eng'] = P2_R3_bond_eng
                break
        ##print("X-R3 bond occ :", P2_R3_bond_occ)
        ##print("X-R3 bond eng :", P2_R3_bond_eng)



    P1_R2_bond_occ = []
    P1_R2_bond_eng = []
    SPE_dict['P1_R2_bond_occ'] = 'Error'
    SPE_dict['P1_R2_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R2'] + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R2']+" - P\s+"+ atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                P1_R2_bond_occ = str.split(line)[9]
                P1_R2_bond_eng = str.split(line)[10]
                SPE_dict['P1_R2_bond_occ'] = P1_R2_bond_occ
                SPE_dict['P1_R2_bond_eng'] = P1_R2_bond_eng
                break
        ##print("P1-R2 bond occ :", P1_R2_bond_occ)
        ##print("P1-R2 bond eng :", P1_R2_bond_eng)

    P2_R4_bond_occ = []
    P2_R4_bond_eng = []
    SPE_dict['P2_R4_bond_occ'] = 'Error'
    SPE_dict['P2_R4_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R4'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R4']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                P2_R4_bond_occ = str.split(line)[9]
                P2_R4_bond_eng = str.split(line)[10]
                SPE_dict['P2_R4_bond_occ'] = P2_R4_bond_occ
                SPE_dict['P2_R4_bond_eng'] = P2_R4_bond_eng
                break
        ##print("X-R4 bond occ :", P2_R4_bond_occ)
        ##print("X-R4 bond eng :", P2_R4_bond_eng)



    P1_RBack1_bond_occ = []
    P1_RBack1_bond_eng = []
    SPE_dict['P1_RBack1_bond_occ'] = 'Error'
    SPE_dict['P1_RBack1_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\)\s+P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['RBack1'] + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack1']+" - P\s+"+ atom_label_dict['P1']
            if re.search(pattern, line):
                P1_RBack1_bond_occ = str.split(line)[9]
                P1_RBack1_bond_eng = str.split(line)[10]
                SPE_dict['P1_RBack1_bond_occ'] = P1_RBack1_bond_occ
                SPE_dict['P1_RBack1_bond_eng'] = P1_RBack1_bond_eng
                break
        ##print("P1-Rback bond occ :", P1_RBack1_bond_occ)
        ##print("P1-Rback bond eng :", P1_RBack1_bond_eng)


    P2_RBack2_bond_occ = []
    P2_Rback_bond_eng = []
    SPE_dict['P2_RBack2_bond_occ'] = 'Error'
    SPE_dict['P2_Rback_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['RBack2'] + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack2']+" - [PN]\s+"+ atom_label_dict['P2']
            if re.search(pattern, line):
                P2_RBack2_bond_occ = str.split(line)[9]
                P2_Rback_bond_eng = str.split(line)[10]
                SPE_dict['P2_RBack2_bond_occ'] = P2_RBack2_bond_occ
                SPE_dict['P2_Rback_bond_eng'] = P2_Rback_bond_eng
                break
        #print("X-Rback bond occ :", P2_RBack2_bond_occ)
        ##print("X-Rback bond eng :", P2_Rback_bond_eng)

    P1_R1_antibond_occ = []
    P1_R1_antibond_eng = []
    SPE_dict['P1_R1_antibond_occ'] = 'Error'
    SPE_dict['P1_R1_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R1'] + '\s+' + '|' + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R1']+" - P\s+"+ atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                P1_R1_antibond_occ = str.split(line)[8]
                P1_R1_antibond_eng = str.split(line)[9]
                SPE_dict['P1_R1_antibond_occ'] = P1_R1_antibond_occ
                SPE_dict['P1_R1_antibond_eng'] = P1_R1_antibond_eng
                break
        ##print("P1-R1 antibond occ :", P1_R1_antibond_occ)
        ##print("P1-R1 antibond eng :", P1_R1_antibond_eng)


    P2_R3_antibond_occ = []
    P2_R3_antibond_eng = []
    SPE_dict['P2_R3_antibond_occ'] = 'Error'
    SPE_dict['P2_R3_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R3'] + '\s+' + '|' + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R3']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                P2_R3_antibond_occ = str.split(line)[8]
                P2_R3_antibond_eng = str.split(line)[9]
                SPE_dict['P2_R3_antibond_occ'] = P2_R3_antibond_occ
                SPE_dict['P2_R3_antibond_eng'] = P2_R3_antibond_eng
                break
        ##print("X-R3 antibond occ :", P2_R3_antibond_occ)
        ##print("X-R3 antibond eng :", P2_R3_antibond_eng)


    P1_R2_antibond_occ = []
    P1_R2_antibond_eng = []
    SPE_dict['P1_R2_antibond_occ'] = 'Error'
    SPE_dict['P1_R2_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R2'] + '\s+' + '|' + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R2']+" - P\s+"+ atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                P1_R2_antibond_occ = str.split(line)[8]
                P1_R2_antibond_eng = str.split(line)[9]
                SPE_dict['P1_R2_antibond_occ'] = P1_R2_antibond_occ
                SPE_dict['P1_R2_antibond_eng'] = P1_R2_antibond_eng
                break
        ##print("P1-R2 antibond occ :", P1_R2_antibond_occ)
        ##print("P1-R2 antibond eng :", P1_R2_antibond_eng)

    P2_R4_antibond_occ = []
    P2_R4_antibond_eng = []
    SPE_dict['P2_R4_antibond_occ'] = 'Error'
    SPE_dict['P2_R4_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R4'] + '\s+' + '|' + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R4']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                P2_R4_antibond_occ = str.split(line)[8]
                P2_R4_antibond_eng = str.split(line)[9]
                SPE_dict['P2_R4_antibond_occ'] = P2_R4_antibond_occ
                SPE_dict['P2_R4_antibond_eng'] = P2_R4_antibond_eng
                break
        ##print("X-R4 antibond occ :", P2_R4_antibond_occ)
        ##print("X-R4 antibond eng :", P2_R4_antibond_eng)



    P1_RBack1_antibond_occ = []
    P1_RBack1_antibond_eng = []
    SPE_dict['P1_RBack1_antibond_occ'] = 'Error'
    SPE_dict['P1_RBack1_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\)\s+P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['RBack1'] + "|" + "BD\*\(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack1']+" - P\s+"+ atom_label_dict['P1']
            if re.search(pattern, line):
                P1_RBack1_antibond_occ = str.split(line)[8]
                P1_RBack1_antibond_eng = str.split(line)[9]
                SPE_dict['P1_RBack1_antibond_occ'] = P1_RBack1_antibond_occ
                SPE_dict['P1_RBack1_antibond_eng'] = P1_RBack1_antibond_eng
                break
        ##print("P1-Rback antibond occ :", P1_RBack1_antibond_occ)
        ##print("P1-Rback antibond eng :", P1_RBack1_antibond_eng)




    P2_RBack2_antibond_occ = []
    P2_RBack2_antibond_eng = []
    SPE_dict['P2_RBack2_antibond_occ'] = 'Error'
    SPE_dict['P2_RBack2_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['RBack2']+ "|" + "BD\*\(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack2']+" - [PN]\s+"+ atom_label_dict['P2']
            if re.search(pattern, line):
                P2_RBack2_antibond_occ = str.split(line)[8]
                P2_RBack2_antibond_eng = str.split(line)[9]
                SPE_dict['P2_RBack2_antibond_occ'] = P2_RBack2_antibond_occ
                SPE_dict['P2_RBack2_antibond_eng'] = P2_RBack2_antibond_eng
                break
        #print("X-Rback antibond occ :", P2_RBack2_antibond_occ)
        #print("X-Rback antibond eng :", P2_RBack2_antibond_eng)



    #----- the below code extracts coordinates from the .log file and returns a dicitonary populated with measurements and angles.

    measurements_dict = {}

    streams, errors = get_outstreams(filename_nolog)


    def get_geom(streams): # extracts the geometry from the compressed stream
        geom = []
        for item in streams[-1][16:]:
            if item == "":
                break
            geom.append([item.split(",")[0],float(item.split(",")[-3]),float(item.split(",")[-2]),float(item.split(",")[-1])])
        return(geom)



    coordinates = get_geom(streams)


    P1_coord = coordinates[int(atom_label_dict['P1'])-1][1:]
    P2_coord = coordinates[int(atom_label_dict['P2'])-1][1:]
    Cl1_coord = coordinates[int(atom_label_dict['Cl1'])-1][1:]
    Cl2_coord = coordinates[int(atom_label_dict['Cl2'])-1][1:]
    Pd_coord = coordinates[int(atom_label_dict['Pd'])-1][1:]
    RBack1_coord = coordinates[int(atom_label_dict['RBack1'])-1][1:]
    R1_coord = coordinates[int(atom_label_dict['R1'])-1][1:]
    R2_coord = coordinates[int(atom_label_dict['R2'])-1][1:]
    RBack2_coord = coordinates[int(atom_label_dict['RBack2'])-1][1:]
    R3_coord = coordinates[int(atom_label_dict['R3'])-1][1:]
    R4_coord = coordinates[int(atom_label_dict['R4'])-1][1:]




    #P1-Pd distance
    x_dist = (P1_coord[0] - Pd_coord[0])**2
    y_dist = (P1_coord[1] - Pd_coord[1])**2
    z_dist = (P1_coord[2] - Pd_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P1-Pd_distance'] = t


    #X_Pd_distance
    x_dist = (P2_coord[0] - Pd_coord[0])**2
    y_dist = (P2_coord[1] - Pd_coord[1])**2
    z_dist = (P2_coord[2] - Pd_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P2-Pd_distance'] = t




    #Cl1_Pd_distance
    x_dist = (Cl1_coord[0] - Pd_coord[0])**2
    y_dist = (Cl1_coord[1] - Pd_coord[1])**2
    z_dist = (Cl1_coord[2] - Pd_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['Cl1-Pd_distance'] = t


    #Cl2_Pd_distance
    x_dist = (Cl2_coord[0] - Pd_coord[0])**2
    y_dist = (Cl2_coord[1] - Pd_coord[1])**2
    z_dist = (Cl2_coord[2] - Pd_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['Cl2-Pd_distance'] = t

    #P1_R1_distance
    x_dist = (P1_coord[0] - R1_coord[0])**2
    y_dist = (P1_coord[1] - R1_coord[1])**2
    z_dist = (P1_coord[2] - R1_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P1_R1_distance'] = t

    #P1_R2_distance
    x_dist = (P1_coord[0] - R2_coord[0])**2
    y_dist = (P1_coord[1] - R2_coord[1])**2
    z_dist = (P1_coord[2] - R2_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P1_R2_distance'] = t

    #P2_R3_distance
    x_dist = (P2_coord[0] - R3_coord[0])**2
    y_dist = (P2_coord[1] - R3_coord[1])**2
    z_dist = (P2_coord[2] - R3_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P2_R3_distance'] = t

    #P2_R4_distance
    x_dist = (P2_coord[0] - R4_coord[0])**2
    y_dist = (P2_coord[1] - R4_coord[1])**2
    z_dist = (P2_coord[2] - R4_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P2_R4_distance'] = t

    #P1_RBack1_distance
    x_dist = (P1_coord[0] - RBack1_coord[0])**2
    y_dist = (P1_coord[1] - RBack1_coord[1])**2
    z_dist = (P1_coord[2] - RBack1_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P1_RBack1_distance'] = t

    #P2_RBack2_distance
    x_dist = (P2_coord[0] - RBack2_coord[0])**2
    y_dist = (P2_coord[1] - RBack2_coord[1])**2
    z_dist = (P2_coord[2] - RBack2_coord[2])**2
    t = math.sqrt(x_dist + y_dist + z_dist)
    measurements_dict['P2_RBack2_distance'] = t

    #bite_angle
    P1 = np.array(P1_coord)
    Pd = np.array(Pd_coord)
    P2 = np.array(P2_coord)

    PdP1 = P1 - Pd
    PdP2 = P2 - Pd
    cosine_angle = np.dot(PdP1, PdP2) / (np.linalg.norm(PdP1)* np.linalg.norm(PdP2))
    angle = np.arccos(cosine_angle)
    measurements_dict['bite_angle'] = np.degrees(angle)



    #distor_P1
    P1 = np.array(P1_coord)
    Pd = np.array(Pd_coord)
    R1 = np.array(R1_coord)
    R2 = np.array(R2_coord)
    RBack1 = np.array(RBack1_coord)

    PdP1 = P1 - Pd
    R1P1 = P1 - R1
    R2P1 = P1 - R2
    RbackP1 = P1 - RBack1

    cosine_angle_R1P1R2 = np.dot(R1P1, R2P1) / (np.linalg.norm(R1P1)* np.linalg.norm(R2P1))
    angle_R1P1R2 = np.arccos(cosine_angle_R1P1R2)

    cosine_angle_R1P1Pd = np.dot(R1P1, PdP1) / (np.linalg.norm(R1P1)* np.linalg.norm(PdP1))
    angle_R1P1Pd = np.arccos(cosine_angle_R1P1Pd)

    cosine_angle_PdP1RBack1 = np.dot(PdP1, RbackP1) / (np.linalg.norm(PdP1)* np.linalg.norm(RbackP1))
    angle_PdP1RBack1 = np.arccos(cosine_angle_PdP1RBack1)

    cosine_angle_R1P1RBack1 = np.dot(R1P1, RbackP1) / (np.linalg.norm(R1P1)* np.linalg.norm(RbackP1))
    angle_R1P1RBack1 = np.arccos(cosine_angle_R1P1RBack1)

    cosine_angle_R2P1RBack1 = np.dot(R2P1, RbackP1) / (np.linalg.norm(R2P1)* np.linalg.norm(RbackP1))
    angle_R2P1RBack1 = np.arccos(cosine_angle_R2P1RBack1)

    cosine_angle_R2P1Pd = np.dot(R2P1, PdP1) / (np.linalg.norm(R2P1)* np.linalg.norm(PdP1))
    angle_R2P1Pd = np.arccos(cosine_angle_R2P1Pd)


    a1= np.degrees(angle_R1P1R2)
    a2= np.degrees(angle_R1P1Pd)
    a3= np.degrees(angle_PdP1RBack1)
    a4= np.degrees(angle_R1P1RBack1)
    a5= np.degrees(angle_R2P1RBack1)
    a6= np.degrees(angle_R2P1Pd)


    measurements_dict['angle_R1P1R2'] =a1
    measurements_dict['angle_R1P1Pd'] =a2
    measurements_dict['angle_PdP1RBack1'] =a3
    measurements_dict['angle_R1P1RBack1'] =a4
    measurements_dict['angle_R2P1RBack1'] =a5
    measurements_dict['angle_R2P1Pd'] =a6

    #distor_X
    P2 = np.array(P2_coord)
    Pd = np.array(Pd_coord)
    R3 = np.array(R3_coord)
    R4 = np.array(R4_coord)
    RBack2 = np.array(RBack2_coord)



    P2Pd = P2 - Pd
    R3P2 = P2 - R3
    R4P2 = P2 - R4
    RBack2 = P2 - RBack2

    cosine_angle_R3P2R4 = np.dot(R3P2, R4P2) / (np.linalg.norm(R3P2)* np.linalg.norm(R4P2))
    angle_R3P2R4 = np.arccos(cosine_angle_R3P2R4)

    cosine_angle_R3P2Pd = np.dot(R3P2, P2Pd) / (np.linalg.norm(R3P2)* np.linalg.norm(P2Pd))
    angle_R3P2Pd = np.arccos(cosine_angle_R3P2Pd)

    cosine_angle_PdP2RBack2 = np.dot(P2Pd, RBack2) / (np.linalg.norm(P2Pd)* np.linalg.norm(RBack2))
    angle_PdP2RBack2 = np.arccos(cosine_angle_PdP2RBack2)

    cosine_angle_R3P2RBack2 = np.dot(R3P2, RBack2) / (np.linalg.norm(R3P2)* np.linalg.norm(RBack2))
    angle_R3P2RBack2 = np.arccos(cosine_angle_R3P2RBack2)

    cosine_angle_R4P2RBack2 = np.dot(R4P2, RBack2) / (np.linalg.norm(R4P2)* np.linalg.norm(RBack2))
    angle_R4P2RBack2 = np.arccos(cosine_angle_R4P2RBack2)

    cosine_angle_R4P2Pd = np.dot(R4P2, P2Pd) / (np.linalg.norm(R4P2)* np.linalg.norm(P2Pd))
    angle_R4P2Pd = np.arccos(cosine_angle_R4P2Pd)


    b1= np.degrees(angle_R3P2R4)
    b2= np.degrees(angle_R3P2Pd)
    b3= np.degrees(angle_PdP2RBack2)
    b4= np.degrees(angle_R3P2RBack2)
    b5= np.degrees(angle_R4P2RBack2)
    b6= np.degrees(angle_R4P2Pd)

    measurements_dict['angle_R3P2R4'] =b1
    measurements_dict['angle_R3P2Pd'] =b2
    measurements_dict['angle_PdP2RBack2'] =b3
    measurements_dict['angle_R3P2RBack2'] =b4
    measurements_dict['angle_R4P2RBack2'] =b5
    measurements_dict['angle_R4P2Pd'] =b6


#Angle_Sum parameters (LvD added in Dec 2021)
    measurements_dict['P1_Angle_Sum'] = a4 + a1 + a5
    measurements_dict['X_Angle_Sum'] = b4 + b1 + b5



    return(atom_label_dict, SPE_dict, measurements_dict)


def get_SPE_NoPd(filename_nolog, atoms_P, atoms_R, atoms_backbone):
    SPE_NoPd_dict = {}
    filename = filename_nolog + '.log'
    atom_keys = ['P1', 'P2', 'R1', 'R2', 'R3', 'R4', 'RBack1', 'RBack2']
    atom_values= atoms_P + atoms_R + atoms_backbone
    atom_label_dict =dict(zip(atom_keys, atom_values))

    HOMO = []
    SPE_NoPd_dict['NoPd_Homo']= 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = r"Alpha\s+occ"
            if re.search(pattern, line):
                HOMO = str.split(line)[-1]
                break
        ##print('NoPd_HOMO:', HOMO)
        SPE_NoPd_dict['NoPd_Homo']= HOMO


    LUMO = []
    SPE_NoPd_dict['NoPd_Lumo'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = r"Alpha\s+virt\."
            if re.search(pattern, line):
                LUMO = str.split(line)[4]
                break
        ##print('LUMO:', LUMO)
        SPE_NoPd_dict['NoPd_Lumo'] = LUMO


    NoPd_dipole = []
    SPE_NoPd_dict['NoPd_dipole'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "\s+X=\s+"
            if re.search(pattern, line):
                NoPd_dipole = str.split(line)[-1]
                break
        ##print("NoPd dipole: ", NoPd_dipole)
        SPE_NoPd_dict['NoPd_dipole'] = NoPd_dipole


    NoPd_R1_NBO = ''
    SPE_NoPd_dict['NoPd_R1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R1']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['R1']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_R1_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_R1_NBO'] = NoPd_R1_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['R1']), filename.split('/')[-1][:-4], NoPd_R1_NBO)

    NoPd_R2_NBO = ''
    SPE_NoPd_dict['NoPd_R2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R2']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['R2']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_R2_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_R2_NBO'] = NoPd_R2_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['R2']), filename.split('/')[-1][:-4], NoPd_R2_NBO)

    NoPd_RBack1_NBO = ''
    SPE_NoPd_dict['NoPd_RBack1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['RBack1']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['RBack1']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_RBack1_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_RBack1_NBO'] = NoPd_RBack1_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['RBack1']), filename.split('/')[-1][:-4], NoPd_RBack1_NBO)

    NoPd_R3_NBO = ''
    SPE_NoPd_dict['NoPd_R3_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R3']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['R3']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_R3_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_R3_NBO'] = NoPd_R3_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['R3']), filename.split('/')[-1][:-4], NoPd_R3_NBO)

    NoPd_R4_NBO = ''
    SPE_NoPd_dict['NoPd_R4_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['R4']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['R4']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_R4_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_R4_NBO'] = NoPd_R4_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['R4']), filename.split('/')[-1][:-4], NoPd_R4_NBO)

    NoPd_RBack2_NBO = ''
    SPE_NoPd_dict['NoPd_RBack2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[CON]\s+" + atom_label_dict['RBack2']+"  \s+-0." + '|' + "[CON]\s+" + atom_label_dict['RBack2']+"  \s+0."
            if re.search(pattern, line):
                if str.split(line)[0] == "C" or str.split(line)[0] == "O" or str.split(line)[0] == "N":
                    NoPd_RBack2_NBO = str.split(line)[2]
                    SPE_NoPd_dict['NoPd_RBack2_NBO'] = NoPd_RBack2_NBO
                break
    #print('anum = {}    '.format(atom_label_dict['RBack2']), filename.split('/')[-1][:-4], NoPd_RBack2_NBO)

    NoPd_P1_NBO = []
    SPE_NoPd_dict['NoPd_P1_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "P\s+" + atom_label_dict['P1']+"    [10]"+ '\.'
            if re.search(pattern, line):
                NoPd_P1_NBO = str.split(line)[2]
                SPE_NoPd_dict['NoPd_P1_NBO'] = NoPd_P1_NBO
                break
        ##print("NoPd_P1 NBO:", NoPd_P1_NBO)


    NoPd_P2_NBO = []
    SPE_NoPd_dict['NoPd_P2_NBO'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "[PN]\s+" + atom_label_dict['P2']+"   [- ][10]"+ '\.'
            if re.search(pattern, line):
                #if str.split(line)[0] == "[PN]":
                NoPd_P2_NBO = str.split(line)[2]
                SPE_NoPd_dict['NoPd_P2_NBO'] = NoPd_P2_NBO
                break
        #print("NoPd_X NBO:", NoPd_P2_NBO)

    NoPd_P1_NMR = ''
    SPE_NoPd_dict['NoPd_P1_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P1']+'\s+P\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    NoPd_P1_NMR = (str.split(line)[4])
                    ##print('NoPd_P1 NMR: ', NoPd_P1_NMR)
                    SPE_NoPd_dict['NoPd_P1_NMR'] = NoPd_P1_NMR
                break


    NoPd_P2_NMR = ''
    SPE_NoPd_dict['NoPd_P2_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P2']+'\s+[PN]\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    NoPd_P2_NMR = (str.split(line)[4])
                    ##print('NoPd_X NMR: ', NoPd_P2_NMR)
                    SPE_NoPd_dict['NoPd_P2_NMR'] = NoPd_P2_NMR
                break


#anisotropic NMR values (added by LvD Nov 2021)
    NoPd_aniso_P1_NMR = ''
    SPE_NoPd_dict['NoPd_aniso_P1_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P1']+'\s+P\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    NoPd_aniso_P1_NMR = (str.split(line)[7])
                    ##print('NoPd_P1 NMR: ', NoPd_P1_NMR)
                    SPE_NoPd_dict['NoPd_aniso_P1_NMR'] = NoPd_aniso_P1_NMR
                break

    NoPd_aniso_P2_NMR = ''
    SPE_NoPd_dict['NoPd_aniso_P2_NMR'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = atom_label_dict['P2']+'\s+[PN]\s+Isotropic'
            if re.search(pattern, line):
                if str.split(line)[2] == 'Isotropic' and str.split(line)[3] == '=':
                    NoPd_aniso_P2_NMR = (str.split(line)[7])
                    #print('X NMR: ', P2_NMR)
                    SPE_NoPd_dict['NoPd_aniso_P2_NMR'] = NoPd_aniso_P2_NMR
                break


    NoPd_P1_LP_occ = []
    NoPd_P1_LP_eng = []
    SPE_NoPd_dict['NoPd_P1_LP_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_LP_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP\s+\(\s+1\)\s+P\s+" + atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                NoPd_P1_LP_occ = str.split(line)[6]
                NoPd_P1_LP_eng = str.split(line)[7]
                SPE_NoPd_dict['NoPd_P1_LP_occ'] = NoPd_P1_LP_occ
                SPE_NoPd_dict['NoPd_P1_LP_eng'] = NoPd_P1_LP_eng
                break
        ##print("NoPd_P1 LP occ :", NoPd_P1_LP_occ)
        ##print("NoPd_P1 LP eng :", NoPd_P1_LP_eng)


    NoPd_P2_LP_occ = []
    NoPd_P2_LP_eng = []
    SPE_NoPd_dict['NoPd_P2_LP_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_LP_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "LP\s+\(\s+1\)\s+[PN]\s+" + atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                NoPd_P2_LP_occ = str.split(line)[6]
                NoPd_P2_LP_eng = str.split(line)[7]
                SPE_NoPd_dict['NoPd_P2_LP_occ'] = NoPd_P2_LP_occ
                SPE_NoPd_dict['NoPd_P2_LP_eng'] = NoPd_P2_LP_eng
                break
        ##print("NoPd_X LP occ :", NoPd_P2_LP_occ)
        ##print("NoPd_X LP eng :", NoPd_P2_LP_eng)


    NoPd_P1_R1_bond_occ = []
    NoPd_P1_R1_bond_eng = []
    SPE_NoPd_dict['NoPd_P1_R1_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_R1_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R1'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R1']+" - P\s+"+ atom_label_dict['P1']+ '\s+'
            if re.search(pattern, line):
                NoPd_P1_R1_bond_occ = str.split(line)[9]
                NoPd_P1_R1_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P1_R1_bond_occ'] = NoPd_P1_R1_bond_occ
                SPE_NoPd_dict['NoPd_P1_R1_bond_eng'] = NoPd_P1_R1_bond_eng
                break
        #print("NoPd_P1-R1 bond occ :", NoPd_P1_R1_bond_occ)
        #print("NoPd_P1-R1 bond eng :", NoPd_P1_R1_bond_eng)


    NoPd_P1_R2_bond_occ = []
    NoPd_P1_R2_bond_eng = []
    SPE_NoPd_dict['NoPd_P1_R2_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_R2_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R2'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R2']+" - P\s+"+ atom_label_dict['P1']+ '\s+'
            if re.search(pattern, line):
                NoPd_P1_R2_bond_occ = str.split(line)[9]
                NoPd_P1_R2_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P1_R2_bond_occ'] = NoPd_P1_R2_bond_occ
                SPE_NoPd_dict['NoPd_P1_R2_bond_eng'] = NoPd_P1_R2_bond_eng
                break
        #print("NoPd_P1-R2 bond occ :", NoPd_P1_R2_bond_occ)
        #print("NoPd_P1-R2 bond eng :", NoPd_P1_R2_bond_eng)


    NoPd_P1_RBack1_bond_occ = []
    NoPd_P1_RBack1_bond_eng = []
    SPE_NoPd_dict['NoPd_P1_RBack1_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_RBack1_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\)\s+P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['RBack1']+ '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack1']+" - P\s+"+ atom_label_dict['P1']+ '\s+'
            if re.search(pattern, line):
                NoPd_P1_RBack1_bond_occ = str.split(line)[9]
                NoPd_P1_RBack1_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P1_RBack1_bond_occ'] = NoPd_P1_RBack1_bond_occ
                SPE_NoPd_dict['NoPd_P1_RBack1_bond_eng'] = NoPd_P1_RBack1_bond_eng
                break
        #print("NoPd_P1-Rback bond occ :", NoPd_P1_RBack1_bond_occ)
        #print("NoPd_P1-Rback bond eng :", NoPd_P1_RBack1_bond_eng)


    NoPd_P2_R3_bond_occ = []
    NoPd_P2_R3_bond_eng = []
    SPE_NoPd_dict['NoPd_P2_R3_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_R3_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R3'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R3']+" - [PN]\s+"+ atom_label_dict['P2']+ '\s+'
            if re.search(pattern, line):
                NoPd_P2_R3_bond_occ = str.split(line)[9]
                NoPd_P2_R3_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P2_R3_bond_occ'] = NoPd_P2_R3_bond_occ
                SPE_NoPd_dict['NoPd_P2_R3_bond_eng'] = NoPd_P2_R3_bond_eng
                break
        #print("NoPd_X_R3 bond occ :", NoPd_P2_R3_bond_occ)
        #print("NoPd_X_R3 bond eng :", NoPd_P2_R3_bond_eng)


    NoPd_P2_R4_bond_occ = []
    NoPd_P2_R4_bond_eng = []
    SPE_NoPd_dict['NoPd_P2_R4_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_R4_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD \(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R4'] + '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['R4']+" - [PN]\s+"+ atom_label_dict['P2']+ '\s+'
            if re.search(pattern, line):
                NoPd_P2_R4_bond_occ = str.split(line)[9]
                NoPd_P2_R4_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P2_R4_bond_occ'] = NoPd_P2_R4_bond_occ
                SPE_NoPd_dict['NoPd_P2_R4_bond_eng'] = NoPd_P2_R4_bond_eng
                break
        #print("NoPd_X_R4 bond occ :", NoPd_P2_R4_bond_occ)
        #print("NoPd_X_R4 bond eng :", NoPd_P2_R4_bond_eng)

    NoPd_P2_RBack2_bond_occ = []
    NoPd_P2_RBack2_bond_eng = []
    SPE_NoPd_dict['NoPd_P2_RBack2_bond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_RBack2_bond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            #pattern = "BD \(\s+1\)\s+P\s+"+ atom_label_dict['P2']+" - C\s+"+ atom_label_dict['RBack2']+ '\s+'
            pattern = "BD \(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['RBack2']+ '\s+' + "|" + "BD \(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack2']+" - [PN]\s+"+ atom_label_dict['P2']+ '\s+'
            if re.search(pattern, line):
                NoPd_P2_RBack2_bond_occ = str.split(line)[9]
                NoPd_P2_RBack2_bond_eng = str.split(line)[10]
                SPE_NoPd_dict['NoPd_P2_RBack2_bond_occ'] = NoPd_P2_RBack2_bond_occ
                SPE_NoPd_dict['NoPd_P2_RBack2_bond_eng'] = NoPd_P2_RBack2_bond_eng
                break
        #print("NoPd_X-Rback bond occ :", NoPd_P2_RBack2_bond_occ)
        #print("NoPd_X-Rback bond eng :", NoPd_P2_RBack2_bond_eng)

##antibonding orbitals
    NoPd_P1_R1_antibond_occ = []
    NoPd_P1_R1_antibond_eng = []
    SPE_NoPd_dict['NoPd_P1_R1_antibond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_R1_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R1'] + '\s+' + "|" + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R1']+" - P\s+"+ atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                NoPd_P1_R1_antibond_occ = str.split(line)[8]
                NoPd_P1_R1_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P1_R1_antibond_occ'] = NoPd_P1_R1_antibond_occ
                SPE_NoPd_dict['NoPd_P1_R1_antibond_eng'] = NoPd_P1_R1_antibond_eng
                break
        #print("NoPd_P1-R1 antibond occ :", NoPd_P1_R1_antibond_occ)
        #print("NoPd_P1-R1 antibond eng :", NoPd_P1_R1_antibond_eng)


    NoPd_P1_R2_antibond_occ = []
    NoPd_P1_R2_antibond_eng = []
    SPE_NoPd_dict['NoPd_P1_R2_antibond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_R2_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['R2'] + '\s+'  +"|"+ "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R2']+" - P\s+"+ atom_label_dict['P1'] + '\s+'
            if re.search(pattern, line):
                NoPd_P1_R2_antibond_occ = str.split(line)[8]
                NoPd_P1_R2_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P1_R2_antibond_occ'] = NoPd_P1_R2_antibond_occ
                SPE_NoPd_dict['NoPd_P1_R2_antibond_eng'] = NoPd_P1_R2_antibond_eng
                break
        #print("NoPd_P1-R2 antibond occ :", NoPd_P1_R2_antibond_occ)
        #print("NoPd_P1-R2 antibond eng :", NoPd_P1_R2_antibond_eng)


    NoPd_P1_RBack1_antibond_occ = []
    NoPd_P1_RBack1_antibond_eng = []
    SPE_NoPd_dict['NoPd_P1_RBack1_antibond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P1_RBack1_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\)\s+P\s+"+ atom_label_dict['P1']+" - [CON]\s+"+ atom_label_dict['RBack1'] + "\s+" +"|"+ "BD\*\(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack1']+" - P\s+"+ atom_label_dict['P1']+ "\s+"
            if re.search(pattern, line):
                NoPd_P1_RBack1_antibond_occ = str.split(line)[8]
                NoPd_P1_RBack1_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P1_RBack1_antibond_occ'] = NoPd_P1_RBack1_antibond_occ
                SPE_NoPd_dict['NoPd_P1_RBack1_antibond_eng'] = NoPd_P1_RBack1_antibond_eng
                break
        #print("NoPd_P1-Rback antibond occ :", NoPd_P1_RBack1_antibond_occ)
        #print("NoPd_P1-Rback antibond eng :", NoPd_P1_RBack1_antibond_eng)


    NoPd_P2_R3_antibond_occ = []
    NoPd_P2_R3_antibond_eng = []
    SPE_NoPd_dict['NoPd_P2_R3_antibond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_R3_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R3'] + '\s+' + "|" + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R3']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                NoPd_P2_R3_antibond_occ = str.split(line)[8]
                NoPd_P2_R3_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P2_R3_antibond_occ'] = NoPd_P2_R3_antibond_occ
                SPE_NoPd_dict['NoPd_P2_R3_antibond_eng'] = NoPd_P2_R3_antibond_eng
                break
        ##print("NoPd_X_R3 antibond occ :", NoPd_P2_R3_antibond_occ)
        ##print("NoPd_X_R3 antibond eng :", NoPd_P2_R3_antibond_eng)

    NoPd_P2_R4_antibond_occ = []
    NoPd_P2_R4_antibond_eng = []
    SPE_NoPd_dict['NoPd_P2_R4_antibond_occ'] = 'Test Error'
    SPE_NoPd_dict['NoPd_P2_R4_antibond_eng'] = 'Test Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\) [PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['R4'] + '\s+' + "|" + "BD\*\(\s+1\) [CON]\s+"+ atom_label_dict['R4']+" - [PN]\s+"+ atom_label_dict['P2'] + '\s+'
            if re.search(pattern, line):
                NoPd_P2_R4_antibond_occ = str.split(line)[8]
                NoPd_P2_R4_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P2_R4_antibond_occ'] = NoPd_P2_R4_antibond_occ
                SPE_NoPd_dict['NoPd_P2_R4_antibond_eng'] = NoPd_P2_R4_antibond_eng
                break
        ##print("NoPd_X_R4 antibond occ :", NoPd_P2_R4_antibond_occ)
        ##print("NoPd_X_R4 antibond eng :", NoPd_P2_R4_antibond_eng)


    NoPd_P2_RBack2_antibond_occ = []
    NoPd_P2_RBack2_antibond_eng = []
    SPE_NoPd_dict['NoPd_P2_RBack2_antibond_occ'] = 'Error'
    SPE_NoPd_dict['NoPd_P2_RBack2_antibond_eng'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f[::-1]:
            pattern = "BD\*\(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+" - [CON]\s+"+ atom_label_dict['RBack2'] + "\s+" +"|"+ "BD\*\(\s+1\)\s+[CON]\s+"+ atom_label_dict['RBack2']+" - [PN]\s+"+ atom_label_dict['P2'] + "\s+"

            if re.search(pattern, line):
                NoPd_P2_RBack2_antibond_occ = str.split(line)[8]
                NoPd_P2_RBack2_antibond_eng = str.split(line)[9]
                SPE_NoPd_dict['NoPd_P2_RBack2_antibond_occ'] = NoPd_P2_RBack2_antibond_occ
                SPE_NoPd_dict['NoPd_P2_RBack2_antibond_eng'] = NoPd_P2_RBack2_antibond_eng
                break
        #print("NoPd_X-Rback antibond occ :", NoPd_P2_RBack2_antibond_occ)
        #print("NoPd_X-Rback antibond eng :", NoPd_P2_RBack2_antibond_eng)


    NoPd_LP_P1_s = []
    SPE_NoPd_dict['NoPd_LP_P1_s'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "\s+LP\s+\(\s+1\)\s+P\s+"+ atom_label_dict['P1']+ "\s+s\("
            if re.search(pattern, line):
                NoPd_LP_P1_s = str.split(line)[8][:-3]
                SPE_NoPd_dict['NoPd_LP_P1_s'] = NoPd_LP_P1_s
                break
        #print("NoPd_LP_P1_s:", NoPd_LP_P1_s)


    NoPd_LP_P2_s = []
    SPE_NoPd_dict['NoPd_LP_P2_s'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "\s+LP\s+\(\s+1\)\s+[PN]\s+"+ atom_label_dict['P2']+ "\s+s\("
            if re.search(pattern, line):
                NoPd_LP_P2_s = str.split(line)[8][:-3]
                SPE_NoPd_dict['NoPd_LP_P2_s'] = NoPd_LP_P2_s
                break
    #print("NoPd_LP_P2_s:", NoPd_LP_P2_s)

    return(SPE_NoPd_dict)

def get_polarizability(filename_nolog):
    polarizability_dict = {}
    filename = filename_nolog + '.log'


    polarizability_dict['pol'] = 'Error'
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "Isotropic polarizability for W="
            if re.search(pattern, line):
                pol = str.split(line)[5]
                ##print("pol:", pol)
                polarizability_dict['pol'] = pol
                break
# Added HF_E and G_energy on  3-10-22 JJD
    HF_E = []
    polarizability_dict['OPT_HF_E'] = 'Error'
    stream = get_outstreams(filename_nolog)
    HF_E = []
    for job in stream[0]:
        for line in job:
            pattern = 'HF=-'
            if re.search(pattern, line):
                HF_E.append(line[3:])
                break
    #print(filename_nolog, HF_E)
    polarizability_dict['OPT_HF_E'] =  HF_E[0]

    polarizability_dict['G'] = 'Error'
    G_energy = []
    with open(filename) as f:
        f = f.readlines()
        for line in f:
            pattern = "thermal Free Energies="
            if re.search(pattern, line):
                G_energy.append(line.split()[-1])
                break
    polarizability_dict['G'] = G_energy[0]
    #print(filename, '  G_energy:  ', G_energy)


    return(polarizability_dict)


def get_outstreams(log): # gets the compressed stream information at the end of a Gaussian job
    streams = []
    starts,ends = [],[]
    error = "failed or incomplete job" # default unless "normal termination" is in file

    try:
        with open(log+".log") as f:
            loglines = f.readlines()
    except:
        with open(log+".LOG") as f:
            loglines = f.readlines()
    for i in range(len(loglines)):
        if "1\\1\\" in loglines[i]:
            starts.append(i)
        if "@" in loglines[i]:
            ends.append(i)
        if "Normal termination" in loglines[i]:
            error = ""
    if len(starts) != len(ends) or len(starts) == 0: #probably redundant
        error = "failed or incomplete job"
        return(streams,error)
    for i in range(len(starts)):
        tmp = ""
        for j in range(starts[i],ends[i]+1,1):
            tmp = tmp + loglines[j][1:-1]
        #print(tmp)
        streams.append(tmp.split("\\"))
    return(streams,error)

def get_geom(streams): # extracts the geometry from the compressed stream
    geom = []
    for item in streams[-1][16:]:
        if item == "":
            break
        geom.append([item.split(",")[0],float(item.split(",")[-3]),float(item.split(",")[-2]),float(item.split(",")[-1])])
    return(geom)

global_bondval_parameters = ['Homo', 'Lumo', 'dipole','Pd_NBO', 'Pd_LP_1_occ', 'Pd_LP_1_eng',
                     'Pd_LP_2_occ', 'Pd_LP_2_eng', 'Pd_LP_3_occ', 'Pd_LP_3_eng',
                     'Pd_LP_4_occ', 'Pd_LP_4_eng', 'bite_angle','NoPd_Homo', 'NoPd_Lumo',
                     'NoPd_dipole','pol']
global_vbur_parameters = ['Vbur%_3.0_Ang','Vbur%_4.0_Ang', 'Vbur%_5.0_Ang', 'Vbur%_6.0_Ang', 'Vbur%_7.0_Ang']
global_parameters = global_bondval_parameters + global_vbur_parameters

C1_dict = {'R1_NBO': ['R1_NBO'], 'R2_NBO': ['R2_NBO'], 'RBack1_NBO': ['RBack1_NBO'], 'R3_NBO': ['R3_NBO'], 'R4_NBO': ['R4_NBO'], 'RBack2_NBO': ['RBack2_NBO'], 'P1_NMR': ['P1_NMR'], 'P2_NMR': ['P2_NMR'], 'aniso_P1_NMR': ['aniso_P1_NMR'], 'aniso_P2_NMR': ['aniso_P2_NMR'], 'P1_NBO': ['P1_NBO'], 'P2_NBO': ['P2_NBO'], 'Cl_1_NBO': ['Cl_1_NBO'], 'Cl_2_NBO': ['Cl_2_NBO'], 'P1_Pd_bond_occ': ['P1_Pd_bond_occ'], 'P2_Pd_bond_occ': ['P2_Pd_bond_occ'], 'P1_Pd_bond_eng': ['P1_Pd_bond_eng'], 'P2_Pd_bond_eng': ['P2_Pd_bond_eng'], 'P1_Pd_antibond_occ': ['P1_Pd_antibond_occ'], 'P2_Pd_antibond_occ': ['P2_Pd_antibond_occ'], 'P1_Pd_antibond_eng': ['P1_Pd_antibond_eng'], 'P2_Pd_antibond_eng': ['P2_Pd_antibond_eng'], 'P1_R1_bond_occ': ['P1_R1_bond_occ'], 'P2_R4_bond_occ': ['P2_R4_bond_occ'], 'P1_R1_bond_eng': ['P1_R1_bond_eng'], 'P2_R4_bond_eng': ['P2_R4_bond_eng'], 'P1_R1_antibond_occ': ['P1_R1_antibond_occ'], 'P2_R4_antibond_occ': ['P2_R4_antibond_occ'], 'P1_R1_antibond_eng': ['P1_R1_antibond_eng'], 'P2_R4_antibond_eng': ['P2_R4_antibond_eng'], 'P1_R2_bond_occ': ['P1_R2_bond_occ'], 'P2_R3_bond_occ': ['P2_R3_bond_occ'], 'P1_R2_bond_eng': ['P1_R2_bond_eng'], 'P2_R3_bond_eng': ['P2_R3_bond_eng'], 'P1_R2_antibond_occ': ['P1_R2_antibond_occ'], 'P2_R3_antibond_occ': ['P2_R3_antibond_occ'], 'P1_R2_antibond_eng': ['P1_R2_antibond_eng'], 'P2_R3_antibond_eng': ['P2_R3_antibond_eng'], 'P1_RBack1_bond_occ': ['P1_RBack1_bond_occ'], 'P2_RBack2_bond_occ': ['P2_RBack2_bond_occ'], 'P1_RBack1_bond_eng': ['P1_RBack1_bond_eng'], 'P2_Rback_bond_eng': ['P2_Rback_bond_eng'], 'P1_RBack1_antibond_occ': ['P1_RBack1_antibond_occ'], 'P2_RBack2_antibond_occ': ['P2_RBack2_antibond_occ'], 'P1_RBack1_antibond_eng': ['P1_RBack1_antibond_eng'], 'P2_RBack2_antibond_eng': ['P2_RBack2_antibond_eng'], 'Cl1-Pd_distance': ['Cl1-Pd_distance'], 'Cl2-Pd_distance': ['Cl2-Pd_distance'], 'P1_R1_distance': ['P1_R1_distance'], 'P1_R2_distance': ['P1_R2_distance'], 'P2_R3_distance': ['P2_R3_distance'], 'P2_R4_distance': ['P2_R4_distance'], 'P1_RBack1_distance': ['P1_RBack1_distance'], 'P2_RBack2_distance': ['P2_RBack2_distance'], 'P1-Pd_distance': ['P1-Pd_distance'], 'P2-Pd_distance': ['P2-Pd_distance'], 'P1_Angle_Sum': ['P1_Angle_Sum'], 'X_Angle_Sum': ['X_Angle_Sum'], 'angle_PdP1RBack1': ['angle_PdP1RBack1'], 'angle_PdP2RBack2': ['angle_PdP2RBack2'], 'angle_R1P1Pd': ['angle_R1P1Pd'], 'angle_R3P2Pd': ['angle_R3P2Pd'], 'angle_R2P1Pd': ['angle_R2P1Pd'], 'angle_R4P2Pd': ['angle_R4P2Pd'], 'angle_R1P1R2': ['angle_R1P1R2'], 'angle_R3P2R4': ['angle_R3P2R4'], 'angle_R1P1RBack1': ['angle_R1P1RBack1'], 'angle_R3P2RBack2': ['angle_R3P2RBack2'], 'angle_R2P1RBack1': ['angle_R2P1RBack1'], 'angle_R4P2RBack2': ['angle_R4P2RBack2'], 'NoPd_R1_NBO': ['NoPd_R1_NBO'], 'NoPd_R2_NBO': ['NoPd_R2_NBO'], 'NoPd_RBack1_NBO': ['NoPd_RBack1_NBO'], 'NoPd_R3_NBO': ['NoPd_R3_NBO'], 'NoPd_R4_NBO': ['NoPd_R4_NBO'], 'NoPd_RBack2_NBO': ['NoPd_RBack2_NBO'], 'NoPd_P1_NMR': ['NoPd_P1_NMR'],
 'NoPd_P2_NMR': ['NoPd_P2_NMR'], 'NoPd_aniso_P1_NMR': ['NoPd_aniso_P1_NMR'], 'NoPd_aniso_P2_NMR': ['NoPd_aniso_P2_NMR'], 'NoPd_P1_NBO': ['NoPd_P1_NBO'], 'NoPd_P2_NBO': ['NoPd_P2_NBO'], 'NoPd_P1_LP_occ': ['NoPd_P1_LP_occ'], 'NoPd_P2_LP_occ': ['NoPd_P2_LP_occ'], 'NoPd_P1_LP_eng': ['NoPd_P1_LP_eng'], 'NoPd_P2_LP_eng': ['NoPd_P2_LP_eng'], 'NoPd_P1_R1_bond_occ': ['NoPd_P1_R1_bond_occ'], 'NoPd_P2_R3_bond_occ': ['NoPd_P2_R3_bond_occ'], 'NoPd_P1_R1_bond_eng': ['NoPd_P1_R1_bond_eng'], 'NoPd_P2_R3_bond_eng': ['NoPd_P2_R3_bond_eng'], 'NoPd_P1_R1_antibond_occ': ['NoPd_P1_R1_antibond_occ'], 'NoPd_P2_R3_antibond_occ': ['NoPd_P2_R3_antibond_occ'], 'NoPd_P1_R1_antibond_eng': ['NoPd_P1_R1_antibond_eng'], 'NoPd_P2_R3_antibond_eng': ['NoPd_P2_R3_antibond_eng'], 'NoPd_P1_R2_bond_occ': ['NoPd_P1_R2_bond_occ'], 'NoPd_P2_R4_bond_occ': ['NoPd_P2_R4_bond_occ'], 'NoPd_P1_R2_bond_eng': ['NoPd_P1_R2_bond_eng'], 'NoPd_P2_R4_bond_eng': ['NoPd_P2_R4_bond_eng'], 'NoPd_P1_R2_antibond_occ': ['NoPd_P1_R2_antibond_occ'], 'NoPd_P2_R4_antibond_occ': ['NoPd_P2_R4_antibond_occ'], 'NoPd_P1_R2_antibond_eng': ['NoPd_P1_R2_antibond_eng'], 'NoPd_P2_R4_antibond_eng': ['NoPd_P2_R4_antibond_eng'], 'NoPd_P1_RBack1_bond_occ': ['NoPd_P1_RBack1_bond_occ'], 'NoPd_P2_RBack2_bond_occ': ['NoPd_P2_RBack2_bond_occ'], 'NoPd_P1_RBack1_bond_eng': ['NoPd_P1_RBack1_bond_eng'], 'NoPd_P2_RBack2_bond_eng': ['NoPd_P2_RBack2_bond_eng'], 'NoPd_P1_RBack1_antibond_occ': ['NoPd_P1_RBack1_antibond_occ'], 'NoPd_P2_RBack2_antibond_occ': ['NoPd_P2_RBack2_antibond_occ'], 'NoPd_P1_RBack1_antibond_eng': ['NoPd_P1_RBack1_antibond_eng'], 'NoPd_P2_RBack2_antibond_eng': ['NoPd_P2_RBack2_antibond_eng'], 'NoPd_LP_P1_s': ['NoPd_LP_P1_s'], 'NoPd_LP_P2_s': ['NoPd_LP_P2_s'], 'NE_3.0_Ang': ['NE_3.0_Ang'], 'SW_3.0_Ang': ['SW_3.0_Ang'], 'NW_3.0_Ang': ['NW_3.0_Ang'], 'SE_3.0_Ang': ['SE_3.0_Ang'], 'NE_4.0_Ang': ['NE_4.0_Ang'], 'SW_4.0_Ang': ['SW_4.0_Ang'], 'NW_4.0_Ang': ['NW_4.0_Ang'], 'SE_4.0_Ang': ['SE_4.0_Ang'], 'NE_5.0_Ang': ['NE_5.0_Ang'], 'SW_5.0_Ang': ['SW_5.0_Ang'], 'NW_5.0_Ang': ['NW_5.0_Ang'], 'SE_5.0_Ang': ['SE_5.0_Ang'],
'SE_6.0_Ang': ['SE_6.0_Ang'],'NE_6.0_Ang': ['NE_6.0_Ang'], 'SW_6.0_Ang': ['SW_6.0_Ang'], 'NW_6.0_Ang': ['NW_6.0_Ang'], 'SE_7.0_Ang': ['SE_7.0_Ang'], 'NE_7.0_Ang': ['NE_7.0_Ang'], 'SW_7.0_Ang': ['SW_7.0_Ang'], 'NW_7.0_Ang': ['NW_7.0_Ang'],
'NE_octant_3.0_Ang': ['NE_octant_3.0_Ang'], 'NW_octant_3.0_Ang': ['NW_octant_3.0_Ang'],
'SE_octant_3.0_Ang': ['SE_octant_3.0_Ang'], 'SW_octant_3.0_Ang': ['SW_octant_3.0_Ang'],
'NE_octant_4.0_Ang': ['NE_octant_4.0_Ang'], 'NW_octant_4.0_Ang': ['NW_octant_4.0_Ang'],
'SE_octant_4.0_Ang': ['SE_octant_4.0_Ang'], 'SW_octant_4.0_Ang': ['SW_octant_4.0_Ang'],
'NE_octant_5.0_Ang': ['NE_octant_5.0_Ang'], 'NW_octant_5.0_Ang': ['NW_octant_5.0_Ang'],
'SE_octant_5.0_Ang': ['SE_octant_5.0_Ang'], 'SW_octant_5.0_Ang': ['SW_octant_5.0_Ang'],
'NE_octant_6.0_Ang': ['NE_octant_6.0_Ang'], 'NW_octant_6.0_Ang': ['NW_octant_6.0_Ang'],
'SE_octant_6.0_Ang': ['SE_octant_6.0_Ang'], 'SW_octant_6.0_Ang': ['SW_octant_6.0_Ang'],
'NE_octant_7.0_Ang': ['NE_octant_7.0_Ang'], 'NW_octant_7.0_Ang': ['NW_octant_7.0_Ang'],
'SE_octant_7.0_Ang': ['SE_octant_7.0_Ang'], 'SW_octant_7.0_Ang': ['SW_octant_7.0_Ang']}


C2_dict_bond_vals = {
# with Pd Parameters
'R1/4_NBO': ['R1_NBO','R4_NBO'],
'R2/3_NBO': ['R2_NBO','R3_NBO'],
'RBack_NBO': ['RBack1_NBO','RBack2_NBO'],
'P_NMR': ['P1_NMR', 'P2_NMR'],
'aniso_P_NMR': ['aniso_P1_NMR', 'aniso_P2_NMR'],
'P_NBO': ['P1_NBO', 'P2_NBO'],
'Cl_NBO': ['Cl_1_NBO', 'Cl_2_NBO'],
'P_Pd_bond_occ': ['P1_Pd_bond_occ', 'P2_Pd_bond_occ'],
'P_Pd_bond_eng': ['P1_Pd_bond_eng', 'P2_Pd_bond_eng'],
'P_Pd_antibond_occ': ['P1_Pd_antibond_occ', 'P2_Pd_antibond_occ'],
'P_Pd_antibond_eng': ['P1_Pd_antibond_eng',  'P2_Pd_antibond_eng'],
'P1_R1/P2_R4_bond_occ': ['P1_R1_bond_occ','P2_R4_bond_occ'],
'P1_R1/P2_R4_bond_eng': ['P1_R1_bond_eng', 'P2_R4_bond_eng'],
'P1_R1/P2_R4_antibond_occ': ['P1_R1_antibond_occ','P2_R4_antibond_occ'],
'P1_R1/P2_R4_antibond_eng': ['P1_R1_antibond_eng', 'P2_R4_antibond_eng'],
'P1_R2/P2_R3_bond_occ': ['P1_R2_bond_occ','P2_R3_bond_occ'],
'P1_R2/P2_R3_bond_eng': ['P1_R2_bond_eng', 'P2_R3_bond_eng'],
'P1_R2/P2_R3_antibond_occ': ['P1_R2_antibond_occ','P2_R3_antibond_occ'],
'P1_R2/P2_R3_antibond_eng': ['P1_R2_antibond_eng', 'P2_R3_antibond_eng'],
'P_RBack_bond_occ': ['P1_RBack1_bond_occ', 'P2_RBack2_bond_occ'],
'P_RBack_bond_eng': ['P1_RBack1_bond_eng',  'P2_Rback_bond_eng'],
'P_RBack_antibond_occ': ['P1_RBack1_antibond_occ', 'P2_RBack2_antibond_occ'],
'P_RBack_antibond_eng': ['P1_RBack1_antibond_eng',  'P2_RBack2_antibond_eng'],
'Cl-Pd_distance': ['Cl1-Pd_distance', 'Cl2-Pd_distance',],
'P1_R1/P2_R4_distance': ['P1_R1_distance', 'P2_R4_distance'],
'P1_R2/P2_R3_distance': ['P1_R2_distance', 'P2_R3_distance'],
'P1_RBack1/P2_RBack2_distance' : ['P1_RBack1_distance', 'P2_RBack2_distance'],

# Geometric Parameters
'P-Pd_distance': ['P1-Pd_distance','P2-Pd_distance'],
'X_Angle_Sum': ['P1_Angle_Sum', 'X_Angle_Sum'],
'angle_Pd_P1_RBack': ['angle_PdP1RBack1', 'angle_PdP2RBack2',],
'angle_R1_P1_Pd/R4_P2_Pd': ['angle_R1P1Pd','angle_R4P2Pd'],
'angle_R2_P1_Pd/R3_P2_Pd': ['angle_R2P1Pd', 'angle_R3P2Pd'],
'angle_R_P_R': ['angle_R1P1R2','angle_R3P2R4',],
'angle_R1_P1_RBack1/R4_P2_RBack2': ['angle_R1P1RBack1', 'angle_R4P2RBack2'],
'angle_R2_P1_RBack1/R3_P2_RBack2': ['angle_R2P1RBack1',  'angle_R3P2RBack2'],

# NoPd Parameters
'NoPd_R1/R4_NBO':['NoPd_R1_NBO','NoPd_R4_NBO',],
'NoPd_R2/3_NBO':['NoPd_R2_NBO','NoPd_R3_NBO'],
'NoPd_RBack_NBO':['NoPd_RBack1_NBO','NoPd_RBack2_NBO'],
'NoPd_P_NMR': ['NoPd_P1_NMR', 'NoPd_P2_NMR'],
'NoPd_aniso_P_NMR': ['NoPd_aniso_P1_NMR', 'NoPd_aniso_P2_NMR'],
'NoPd_P_NBO': ['NoPd_P1_NBO', 'NoPd_P2_NBO'],
'NoPd_P_LP_occ': ['NoPd_P1_LP_occ','NoPd_P2_LP_occ',],
'NoPd_P_LP_eng': ['NoPd_P1_LP_eng',  'NoPd_P2_LP_eng',],
'NoPd_P1_R1/P2_R4_bond_occ': ['NoPd_P1_R1_bond_occ','NoPd_P2_R4_bond_occ'],
'NoPd_P1_R1/P2_R4_bond_eng': ['NoPd_P1_R1_bond_eng','NoPd_P2_R4_bond_eng'],
'NoPd_P1_R1/P2_R4_antibond_occ': ['NoPd_P1_R1_antibond_occ','NoPd_P2_R4_antibond_occ'],
'NoPd_P1_R1/P2_R4_antibond_eng': ['NoPd_P1_R1_antibond_eng','NoPd_P2_R4_antibond_eng'],
'NoPd_P1_R2/P2_R3_bond_occ': ['NoPd_P1_R2_bond_occ', 'NoPd_P2_R3_bond_occ'],
'NoPd_P1_R2/P2_R3_bond_eng': ['NoPd_P1_R2_bond_eng','NoPd_P2_R3_bond_eng'],
'NoPd_P1_R2/P2_R3_antibond_occ': ['NoPd_P1_R2_antibond_occ', 'NoPd_P2_R3_antibond_occ'],
'NoPd_P1_R2/P2_R3_antibond_eng': ['NoPd_P1_R2_antibond_eng', 'NoPd_P2_R3_antibond_eng'],
'NoPd_P_RBack_bond_occ': ['NoPd_P1_RBack1_bond_occ','NoPd_P2_RBack2_bond_occ'],
'NoPd_P_RBack_bond_eng': ['NoPd_P1_RBack1_bond_eng', 'NoPd_P2_RBack2_bond_eng'],
'NoPd_P_RBack_antibond_occ': ['NoPd_P1_RBack1_antibond_occ', 'NoPd_P2_RBack2_antibond_occ'],
'NoPd_P_RBack_antibond_eng': ['NoPd_P1_RBack1_antibond_eng', 'NoPd_P2_RBack2_antibond_eng'],
'NoPd_LP_P_s': ['NoPd_LP_P1_s', 'NoPd_LP_P2_s']}


C2_dict_vbur = { 'NE_SW_quadrants_3.0': ['NE_3.0_Ang', 'SW_3.0_Ang',],
                 'NW_SE_quadrants_3.0': ['NW_3.0_Ang', 'SE_3.0_Ang',],
                 'NE_SW_quadrants_4.0': ['NE_4.0_Ang', 'SW_4.0_Ang',],
                 'NW_SE_quadrants_4.0': ['NW_4.0_Ang', 'SE_4.0_Ang',],
                 'NE_SW_quadrants_5.0': ['NE_5.0_Ang', 'SW_5.0_Ang',],
                 'NW_SE_quadrants_5.0': ['NW_5.0_Ang', 'SE_5.0_Ang'],
                 'NE_SW_quadrants_6.0': ['NE_6.0_Ang', 'SW_6.0_Ang',],
                 'NW_SE_quadrants_6.0': ['NW_6.0_Ang', 'SE_6.0_Ang'],
                 'NE_SW_quadrants_7.0': ['NE_7.0_Ang', 'SW_7.0_Ang',],
                 'NW_SE_quadrants_7.0': ['NW_7.0_Ang', 'SE_7.0_Ang'],
                 'NE_SW_octants_3.0_Ang': ['NE_octant_3.0_Ang','SW_octant_3.0_Ang'],
                 'SE_NW_octants_3.0_Ang': ['SE_octant_3.0_Ang','NW_octant_3.0_Ang'],
                 'NE_SW_octants_4.0_Ang': ['NE_octant_4.0_Ang','SW_octant_4.0_Ang'],
                 'SE_NW_octants_4.0_Ang': ['SE_octant_4.0_Ang', 'NW_octant_4.0_Ang'],
                 'NE_SW_octants_5.0_Ang': ['NE_octant_5.0_Ang','SW_octant_5.0_Ang'],
                 'SE_NW_octants_5.0_Ang': ['SE_octant_5.0_Ang','NW_octant_5.0_Ang'],
                 'NE_SW_octants_6.0_Ang': ['NE_octant_6.0_Ang','SW_octant_6.0_Ang'],
                 'SE_NW_octants_6.0_Ang': ['SE_octant_6.0_Ang','NW_octant_6.0_Ang'],
                 'NE_SW_octants_7.0_Ang': ['NE_octant_7.0_Ang','SW_octant_7.0_Ang'],
                 'SE_NW_octants_7.0_Ang': ['SE_octant_7.0_Ang','NW_octant_7.0_Ang']}

C2_dict = {**C2_dict_bond_vals, **C2_dict_vbur}


C2v_dict_bond_vals = {
# With Pd Parameters
'R_NBO': ['R1_NBO','R2_NBO','R3_NBO','R4_NBO'],
'RBack_NBO': ['RBack1_NBO','RBack2_NBO'],
'P_NMR': ['P1_NMR', 'P2_NMR'],
'aniso_P_NMR': ['aniso_P1_NMR', 'aniso_P2_NMR'],
'P_NBO': ['P1_NBO', 'P2_NBO'],
'Cl_NBO': ['Cl_1_NBO', 'Cl_2_NBO'],
'P_Pd_bond_occ': ['P1_Pd_bond_occ', 'P2_Pd_bond_occ'],
'P_Pd_bond_eng': ['P1_Pd_bond_eng', 'P2_Pd_bond_eng'],
'P_Pd_antibond_occ': ['P1_Pd_antibond_occ', 'P2_Pd_antibond_occ'],
'P_Pd_antibond_eng': ['P1_Pd_antibond_eng',  'P2_Pd_antibond_eng'],
'P_R_bond_occ': ['P1_R1_bond_occ','P2_R3_bond_occ', 'P1_R2_bond_occ','P2_R4_bond_occ'],
'P_R_bond_eng': ['P1_R1_bond_eng',  'P2_R3_bond_eng', 'P1_R2_bond_eng',  'P2_R4_bond_eng'],
'P_R_antibond_occ': ['P1_R1_antibond_occ','P2_R3_antibond_occ','P1_R2_antibond_occ','P2_R4_antibond_occ'],
'P_R_antibond_eng': ['P1_R1_antibond_eng',  'P2_R3_antibond_eng','P1_R2_antibond_eng',  'P2_R4_antibond_eng'],
'P_RBack_bond_occ': ['P1_RBack1_bond_occ', 'P2_RBack2_bond_occ'],
'P_RBack_bond_eng': ['P1_RBack1_bond_eng',  'P2_Rback_bond_eng'],
'P_RBack_antibond_occ': ['P1_RBack1_antibond_occ', 'P2_RBack2_antibond_occ'],
'P_RBack_antibond_eng': ['P1_RBack1_antibond_eng',  'P2_RBack2_antibond_eng'],
'P_R_distance': ['P1_R1_distance', 'P2_R4_distance','P1_R2_distance', 'P2_R3_distance'],
'P1_RBack1/P2_RBack2_distance' : ['P1_RBack1_distance', 'P2_RBack2_distance'],

# Geometric Parameters
'Cl-Pd_distance': ['Cl1-Pd_distance', 'Cl2-Pd_distance',],
'P-Pd_distance': ['P1-Pd_distance','P2-Pd_distance'],
'X_Angle_Sum': ['P1_Angle_Sum', 'X_Angle_Sum'],
'angle_Pd_P1_RBack1/Pd_P2_RBack2': ['angle_PdP1RBack1', 'angle_PdP2RBack2',],
'angle_RPPd': ['angle_R1P1Pd','angle_R3P2Pd', 'angle_R2P1Pd', 'angle_R4P2Pd'],
'angle_R1_P1_R2/R3_P2_R4': ['angle_R1P1R2','angle_R3P2R4',],
'angle_RPRBack': ['angle_R1P1RBack1',   'angle_R3P2RBack2', 'angle_R2P1RBack1',  'angle_R4P2RBack2'],

# NoPd Parameters
'NoPd_R_NBO' : ['NoPd_R1_NBO','NoPd_R2_NBO','NoPd_R3_NBO','NoPd_R4_NBO'],
'NoPd_Rback_NBO':['NoPd_RBack1_NBO','NoPd_RBack2_NBO'],
'NoPd_P_NMR': ['NoPd_P1_NMR', 'NoPd_P2_NMR'],
'NoPd_aniso_P_NMR': ['NoPd_aniso_P1_NMR', 'NoPd_aniso_P2_NMR'],
'NoPd_P_NBO': ['NoPd_P1_NBO', 'NoPd_P2_NBO'],
'NoPd_P_LP_occ': ['NoPd_P1_LP_occ','NoPd_P2_LP_occ',],
'NoPd_P_LP_eng': ['NoPd_P1_LP_eng',  'NoPd_P2_LP_eng',],
'NoPd_P_R_bond_occ': ['NoPd_P1_R1_bond_occ','NoPd_P2_R3_bond_occ', 'NoPd_P1_R2_bond_occ', 'NoPd_P2_R4_bond_occ'],
'NoPd_P_R_bond_eng': ['NoPd_P1_R1_bond_eng','NoPd_P2_R3_bond_eng', 'NoPd_P1_R2_bond_eng','NoPd_P2_R4_bond_eng'],
'NoPd_P_R_antibond_occ': ['NoPd_P1_R1_antibond_occ','NoPd_P2_R3_antibond_occ', 'NoPd_P1_R2_antibond_occ', 'NoPd_P2_R4_antibond_occ'],
'NoPd_P_R_antibond_eng': ['NoPd_P1_R1_antibond_eng','NoPd_P2_R3_antibond_eng', 'NoPd_P1_R2_antibond_eng', 'NoPd_P2_R4_antibond_eng'],
'NoPd_P_RBack_bond_occ': ['NoPd_P1_RBack1_bond_occ','NoPd_P2_RBack2_bond_occ'],
'NoPd_P_RBack_bond_eng': ['NoPd_P1_RBack1_bond_eng', 'NoPd_P2_RBack2_bond_eng'],
'NoPd_P_RBack_antibond_occ': ['NoPd_P1_RBack1_antibond_occ', 'NoPd_P2_RBack2_antibond_occ'],
'NoPd_P_RBack_antibond_eng': ['NoPd_P1_RBack1_antibond_eng', 'NoPd_P2_RBack2_antibond_eng'],
'NoPd_LP_P_s': ['NoPd_LP_P1_s', 'NoPd_LP_P2_s']}


C2v_dict_vbur = {'quadrants_3.0_Ang': ['NE_3.0_Ang', 'SW_3.0_Ang','NW_3.0_Ang', 'SE_3.0_Ang'],
                 'quadrants_4.0_Ang': ['NE_4.0_Ang', 'SW_4.0_Ang','NW_4.0_Ang', 'SE_4.0_Ang'],
                 'quadrants_5.0_Ang': ['NE_5.0_Ang', 'SW_5.0_Ang','NW_5.0_Ang', 'SE_5.0_Ang'],
                 'quadrants_6.0_Ang': ['NE_6.0_Ang', 'SW_6.0_Ang','NW_6.0_Ang', 'SE_6.0_Ang'],
                 'quadrants_7.0_Ang': ['NE_7.0_Ang', 'SW_7.0_Ang','NW_7.0_Ang', 'SE_7.0_Ang'],
                 'octants_3.0_Ang': ['NE_octant_3.0_Ang','SW_octant_3.0_Ang', 'SE_octant_3.0_Ang','NW_octant_3.0_Ang'],
                 'octants_4.0_Ang': ['NE_octant_4.0_Ang','SW_octant_4.0_Ang', 'SE_octant_4.0_Ang', 'NW_octant_4.0_Ang'],
                 'octants_5.0_Ang': ['NE_octant_5.0_Ang','SW_octant_5.0_Ang','SE_octant_5.0_Ang','NW_octant_5.0_Ang'],
                 'octants_6.0_Ang': ['NE_octant_6.0_Ang','SW_octant_6.0_Ang','SE_octant_6.0_Ang','NW_octant_6.0_Ang'],
                 'octants_7.0_Ang': ['NE_octant_7.0_Ang','SW_octant_7.0_Ang','SE_octant_7.0_Ang','NW_octant_7.0_Ang']}

C2v_dict = {**C2v_dict_bond_vals, **C2v_dict_vbur}

ratio_dict = {'RATIO_angle_R1-X-Pd/R2-X-Pd': [['angle_R1P1Pd','angle_R2P1Pd'], ['angle_R4P2Pd','angle_R3P2Pd']],
'RATIO_angle_R1-X-Rback/R2-X-Rback': [['angle_R1P1RBack1','angle_R2P1RBack1'], ['angle_R4P2RBack2','angle_R3P2RBack2']],
'RATIO_N/S_quadrant_3.0_Ang': [['NW_3.0_Ang','SW_3.0_Ang'], ['SE_3.0_Ang','NE_3.0_Ang']],
'RATIO_N/S_quadrant_4.0_Ang': [['NW_4.0_Ang','SW_4.0_Ang'], ['SE_4.0_Ang','NE_4.0_Ang']],
'RATIO_N/S_quadrant_5.0_Ang': [['NW_5.0_Ang','SW_5.0_Ang'], ['SE_5.0_Ang','NE_5.0_Ang']],
'RATIO_N/S_octant_3.0_Ang': [['NW_octant_3.0_Ang','SW_octant_3.0_Ang'], ['SE_octant_3.0_Ang','NE_octant_3.0_Ang']],
'RATIO_N/S_octant_4.0_Ang': [['NW_octant_4.0_Ang','SW_octant_4.0_Ang'], ['SE_octant_4.0_Ang','NE_octant_4.0_Ang']],
'RATIO_N/S_octant_5.0_Ang': [['NW_octant_5.0_Ang','SW_octant_5.0_Ang'], ['SE_octant_5.0_Ang','NE_octant_5.0_Ang']],
'RATIO_N/S_octant_6.0_Ang': [['NW_octant_6.0_Ang','SW_octant_6.0_Ang'], ['SE_octant_6.0_Ang','NE_octant_6.0_Ang']],
'RATIO_N/S_octant_7.0_Ang': [['NW_octant_7.0_Ang','SW_octant_7.0_Ang'], ['SE_octant_7.0_Ang','NE_octant_7.0_Ang']]
}
