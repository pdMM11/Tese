# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 15:51:46 2020

@author: pedro
"""
count_total = 0
connt_diff = 0
peri = 0
for i in results_output.keys():
    for j in results_output[i].keys():
        if results_output[i][j] != results_output_fast[i][j]:
            connt_diff += 1
            print(i,j)
            print(results_output[i][j] - results_output_fast[i][j])
            print()
            if j== 'Peribunyaviridae':
                peri += 1
        count_total += 1
print("Error:", connt_diff/count_total)
print("Errors in Peribunyaviridae:", peri / connt_diff) 

"""
arenaviridae41 Peribunyaviridae
0.047512872659370486

arenaviridae141 Peribunyaviridae
0.03851739130434778

bunyaviridae82 Herpesviridae
0.007062700714904391

coronaviridae184 Peribunyaviridae
0.022482477320087624

filoviridae217 Peribunyaviridae
0.015920260378305295

flaviviridae73 Peribunyaviridae
0.042110007413075135

flaviviridae76 Peribunyaviridae
0.014617503706537582

flaviviridae135 Coronaviridae
0.021622704607615062

flaviviridae214 Coronaviridae
0.015437047986961328

flaviviridae212 Coronaviridae
0.021622704607615062

flaviviridae210 Coronaviridae
0.021622704607615062

flaviviridae192 Coronaviridae
0.021622704607615062

flaviviridae163 Coronaviridae
0.021622704607615062

flaviviridae156 Coronaviridae
0.021622704607615062

flaviviridae54 Peribunyaviridae
0.06688186222610876

flaviviridae29 Peribunyaviridae
0.14929687069819464

flaviviridae28 Peribunyaviridae
0.04149487351738501

herpesviridae128 Peribunyaviridae
0.041466404744368135

orthomyxoviridae117 Pneumoviridae
0.01055631957827835

orthomyxoviridae207 Pneumoviridae
0.015073083601001291

orthomyxoviridae123 Pneumoviridae
0.0003735280454457146

orthomyxoviridae157 Pneumoviridae
0.009316653982283207

orthomyxoviridae25 Peribunyaviridae
0.12944193328423076

paramyxoviridae130 Peribunyaviridae
0.044310561897747236

paramyxoviridae140 Peribunyaviridae
0.012164000000000008

paramyxoviridae137 Peribunyaviridae
0.002933280948873629

paramyxoviridae66 Peribunyaviridae
0.0011672809488735836

paramyxoviridae14 Peribunyaviridae
0.030767842846620874

paramyxoviridae218 Peribunyaviridae
0.03128009636770673

paramyxoviridae196 Peribunyaviridae
0.015284671185112508

paramyxoviridae179 Peribunyaviridae
0.050925280948873664

paramyxoviridae5 Peribunyaviridae
0.03416384284662094

paramyxoviridae47 Peribunyaviridae
0.009548000000000001

paramyxoviridae101 Peribunyaviridae
0.06981056189774726

poxviridae63 Peribunyaviridae
0.014984297627815946

retroviridae16 Peribunyaviridae
0.027790102015763984

retroviridae188 Peribunyaviridae
0.00879761904761911

retroviridae96 Peribunyaviridae
0.00960133785177919

retroviridae161 Peribunyaviridae
0.015074478449699169

retroviridae145 Coronaviridae
0.00767580924140665

retroviridae35 Peribunyaviridae
0.030694194994636326

rhabdoviridae168 Peribunyaviridae
0.06266425593046021

rhabdoviridae202 Peribunyaviridae
0.06684187299249089

togaviridae91 Peribunyaviridae
0.06502777777777774

togaviridae99 Peribunyaviridae
0.0016212235401021902

togaviridae19 Peribunyaviridae
0.03580989416040897

togaviridae33 Peribunyaviridae
0.06675672617586237

Error: 0.016285516285516284
Errors in Peribunyaviridae: 0.723404255319149
"""


count_total = 0
connt_diff = 0
peri = 0
for i in results_output.keys():
    for j in results_output[i].keys():
        if j == 'Peribunyaviridae': pass
        else:
            if results_output[i][j] != results_output_fast[i][j]:
                connt_diff += 1
                print(i,j)
                print(results_output[i][j] - results_output_fast[i][j])
                print()
            count_total += 1
print("Error:", connt_diff/count_total)
print("Errors in Peribunyaviridae: ", peri / connt_diff)


"""
bunyaviridae82 Herpesviridae
0.007062700714904391

flaviviridae135 Coronaviridae
0.021622704607615062

flaviviridae214 Coronaviridae
0.015437047986961328

flaviviridae212 Coronaviridae
0.021622704607615062

flaviviridae210 Coronaviridae
0.021622704607615062

flaviviridae192 Coronaviridae
0.021622704607615062

flaviviridae163 Coronaviridae
0.021622704607615062

flaviviridae156 Coronaviridae
0.021622704607615062

orthomyxoviridae117 Pneumoviridae
0.01055631957827835

orthomyxoviridae207 Pneumoviridae
0.015073083601001291

orthomyxoviridae123 Pneumoviridae
0.0003735280454457146

orthomyxoviridae157 Pneumoviridae
0.009316653982283207

retroviridae145 Coronaviridae
0.00767580924140665

Error: 0.00487987987987988
Errors in Peribunyaviridae:  0.0
"""