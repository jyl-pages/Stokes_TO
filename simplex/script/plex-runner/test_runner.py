import os
import sys
import pickle

def param_command(param_dict): # param_dict is like {'-s':80, '-lf':300...}
    L=list(zip(param_dict.keys(),param_dict.values()))
    L1=list(map(lambda tp:str(tp[0])+' '+str(tp[1]),L))
    L2=' '.join(L1)
    return L2

class Test:
    def __init__(self,exam_exe,test_folder,param_dict):
        self.exam_exe=exam_exe
        param_dict['-o']=test_folder
        self.command=param_command(param_dict)
        return

    def run_test(self):
        run_cmd="{} {}".format(self.exam_exe,self.command)
        print(run_cmd)
        os.system(run_cmd)
        return

if __name__=="__main__":
    run_file=sys.argv[1]
    with open(run_file,'rb') as f:
        tests=pickle.load(f)
    N=len(tests)
    bn,k=map(int,sys.argv[2:4])
    for i in range(N):
        #k/bn <= i/N < (k+1)/bn
        #k*N <= i*bn <(k+1)*N
        if k*N <= i*bn and i*bn < (k+1)*N:
            tests[i].run_test()
