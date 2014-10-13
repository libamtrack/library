import logging
import argparse

import src.common as common
import src.generation_library as generation_library 
import src.in_out as in_out
import src.plotting as plotting
import src.pyamtrack
import src.pyamtrack_SPC

from src.in_out import *
from src.generation_library import *
from src.plotting import *

import numpy as np

def printInitSpectrum(positions, coefficients, maximum, distal):
        # others 
    material_no = 1
    stopping_power_source_no = 0

    # scan SPC headers and read pairs of energy + maximum position
    dirname = '../spc/'    
    E_v = []
    p_v = []
    for filename in os.listdir(dirname):
        if filename.endswith('.spc'):
            #print filename
            E_MeV_u = 0.
            peak_position_g_cm2 = 0.
            particle_no = 0
            material_no = 0
            normalisation = 1.
            steps_no = 1
    
            result = src.pyamtrack.AT_SPC_read_header_from_filename_fast(dirname + filename, E_MeV_u, peak_position_g_cm2, particle_no, material_no, normalisation, steps_no)
            E_v.append( result[1] )
            p_v.append( result[2] )
 
    E_v, p_v = zip(*sorted(zip(E_v, p_v)))
    # remove point corresponding to 500MeV , wrong one ?
    E_v = E_v[0:-1]
    p_v = p_v[0:-1]

    # read spectrum corresponding to maximum distal BP    
    spectrum_full = src.pyamtrack_SPC.extract_spectrum_at_range(dirname, distal)
    spectrum = spectrum_full
    spectrum_dict = src.pyamtrack_SPC.get_spectrum_dictionary_depth(spectrum)
    
    # rescale to maximum dose of single BP = 1
    maximum_dose_Gy = 1.
    scaled_spectrum_dict = src.pyamtrack_SPC.scale_fluence_to_maximum_dose( maximum_dose_Gy, maximum, spectrum_dict)
    
    #input SOBP dose
    input_dose_cm2 = sum([coefficients[i] * src.pyamtrack_SPC.dose_at_depth(0, maximum, positions[i], scaled_spectrum_dict) for i in range(len(coefficients)) ])
    #print "input dose = ", input_dose_cm2
    
    #input SOBP fluence        
    input_fluence_cm2 = sum([coefficients[i] * src.pyamtrack_SPC.fluence_at_depth(0, maximum, positions[i], scaled_spectrum_dict)[6012] for i in range(len(coefficients)) ])
    #print "input fluence = ", input_fluence_cm2
    
    components_fluence = [src.pyamtrack_SPC.fluence_at_depth(0, maximum, positions[i], scaled_spectrum_dict)[6012] for i in range(len(coefficients)) ]
    components_dose = [src.pyamtrack_SPC.dose_at_depth(0, maximum, positions[i], scaled_spectrum_dict) for i in range(len(coefficients)) ]
    #print "components fluence = ", components_fluence
    #print "components dose = ", components_dose
    
    for i in range(len(positions)):
        E = np.interp( positions[i], p_v , E_v)
        print positions[i], coefficients[i], E, coefficients[i] * components_fluence[i]

    
    pass

if __name__ == "__main__":

    # positions and coefficients read from generation file
    positions = [8.0, 8.081632653061224, 8.16326530612245, 8.244897959183675, 8.326530612244898, 8.408163265306122, 8.489795918367347, 8.571428571428571, 8.653061224489797, 8.73469387755102, 8.816326530612246, 8.89795918367347, 8.979591836734695, 9.061224489795919, 9.142857142857142, 9.224489795918368, 9.306122448979592, 9.387755102040817, 9.46938775510204, 9.551020408163266, 9.63265306122449, 9.714285714285715, 9.795918367346939, 9.877551020408163, 9.959183673469388, 10.040816326530612, 10.122448979591837, 10.204081632653061, 10.285714285714286, 10.36734693877551, 10.448979591836736, 10.53061224489796, 10.612244897959183, 10.693877551020408, 10.775510204081632, 10.857142857142858, 10.938775510204081, 11.020408163265307, 11.10204081632653, 11.183673469387756, 11.26530612244898, 11.346938775510203, 11.428571428571429, 11.510204081632653, 11.591836734693878, 11.673469387755102, 11.755102040816327, 11.83673469387755, 11.918367346938776, 12.0]
    coefficients =  [ 0.0379807314799 , 0.0324058208475 , 0.034634791914 , 0.0333481098966 , 0.0340855588217 , 0.0348233244977 , 0.0349203392919 , 0.0345865975699 , 0.0345357699411 , 0.0353402992539 , 0.0359827625809 , 0.037110642906 , 0.0379562484903 , 0.0389049003913 , 0.0393642645769 , 0.0395485284457 , 0.0401753267803 , 0.0405154100425 , 0.0413774201366 , 0.0419946213427 , 0.0434639688209 , 0.0442232724525 , 0.0452018386211 , 0.0465301609255 , 0.0477490067816 , 0.049187285352 , 0.0503454022971 , 0.0521237966446 , 0.0532895071478 , 0.0547639671183 , 0.055776157223 , 0.0571825507628 , 0.0585631503387 , 0.0598787995976 , 0.0624565738173 , 0.0650951081967 , 0.0689695969544 , 0.0735744219735 , 0.0781133898366 , 0.0814659637587 , 0.0849737446943 , 0.0882866070695 , 0.0956859830155 , 0.10997453011 , 0.131733165931 , 0.145281550548 , 0.135018770742 , 0.242996386667 , 0.23577357119 , 0.802739150869 ,  ]     
    # max Bragg peak positions
    maximum = 11.77701
    distal = 12
    
    printInitSpectrum( positions, coefficients , maximum, distal)
    
