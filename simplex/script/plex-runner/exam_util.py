import os
import sys
import numpy as np

origin_root=os.path.join('..','..','..')
simplex_root=os.path.join(origin_root,'simplex')
complex_root=os.path.join(origin_root,'complex')

#plex_name: "simplex"/"complex", distr_name: "Debug"/"Release"
def plex_bin_path(plex_name,proj_name,distr_name='Release'):
    plex_root=os.path.join(origin_root,plex_name)
    bin_root=os.path.join(plex_root,'bin','win',proj_name,distr_name)
    bin_path=os.path.join(bin_root,'{}.exe'.format(proj_name))
    return bin_path

def copy_file(file_a,file_b):
    copy_cmd="copy {} {}".format(file_a,file_b)
    os.system(copy_cmd)
    return

def read_int(counter_txt):
    with open(counter_txt,"rt") as fin:
        str=fin.read()
        num=int(str)
    return num

def write_int(counter_txt,num):
    with open(counter_txt,"wt") as fout:
        fout.write(str(num))
    return

class Exam:
    def __init__(self,plex_name,proj_name,distr_name,data_folder,exam_id):
        self.compile_exe=plex_bin_path(plex_name,proj_name,distr_name)
        self.data_folder=data_folder
        self.exam_id=exam_id
        self.test_counter=0
        self.prep_exam_env()
        return

    def prep_exam_env(self):
        self.exam_folder=os.path.join(self.data_folder,self.exam_id)
        os.system("mkdir {}".format(self.exam_folder))
        self.exam_exe=os.path.join(self.exam_folder,"{}.exe".format(self.exam_id))
        copy_file(self.compile_exe,self.exam_exe)
        return

    def prep_test_env(self,test_id):
        test_seq=self.test_counter
        self.test_counter+=1
        test_name="{}-{}".format(test_seq,test_id)
        test_folder=os.path.join(self.exam_folder,test_name)
        os.system("mkdir {}".format(test_folder))
        return test_folder
