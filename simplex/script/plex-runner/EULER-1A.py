import sys
import os
import pickle
import exam_util
import test_runner

if __name__=='__main__':
    #The name of exam is taken as the name of this file, i.e., 'EULER-1A'
    exam_id=__file__.split('.')[0]
    #Specify in script a file that stores required information for running all tests of this exam
    run_file=sys.argv[1]
    #The folder that stores data of all tests
    data_folder=os.path.join('.','data')
    #The project of this exam is simplex/fluid_euler
    #You can change it to any project in simplex or complex freely
    E=exam_util.Exam('simplex','fluid_euler','Release',data_folder,exam_id)
    #Common parameters. '-o' is specially not required here
    param_dict={'-driver':1,'-test':1,'-lf':200}

    tests = []
    for s in [50,100,150,200,250]:
        test_id="s_{}".format(s)
        test_folder=E.prep_test_env(test_id)
        #This exam includes 3 tests with different resolution
        param_dict['-s']=s
        t=test_runner.Test(E.exam_exe,test_folder,param_dict)
        tests.append(t)

    with open(run_file,'wb') as f:
        pickle.dump(tests,f)
