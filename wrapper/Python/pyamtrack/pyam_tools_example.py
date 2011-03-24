import pyam_tools as pyam

# create PyAmRun object, called TEST
test = pyam.PyAmRun('TEST')

# define Attributs of the object
test.algorithm = 'IGK'
test.number_of_components = 1
test.E_MeV_u = [10.0]
test.particle_no = [1001]
test.fluence_cm2 = [-10.0]
test.material_no = 5
test.rdd_model = 4
test.rdd_parameters[:4]= [1e-10,1e-10,0.0, 0.0]
test.er_model =2
test.gamma_model = 2
test.gamma_parameters[:4]=[1,10.5e4,1,1]


# Run algorithm on object
test.run()
igk_re = test.relative_efficiency

# change algorithm
test.algorithm= 'CPPSC'
test.run()
cppsc_re = test.relative_efficiency

# print the results
print '### TEST ###'
print 'RE\nIGK :\t%.3f\nCPPSC :\t%.3f'%(igk_re, cppsc_re)

