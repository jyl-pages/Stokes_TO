from shutil import copyfile
import os
import sys, glob, time

if __name__ == "__main__":
	input_folder = sys.argv[1]; output_folder = sys.argv[2];
	target_file = sys.argv[3]; last_frame = int(sys.argv[4]);
	dst_folder = output_folder;
	if not os.path.exists(dst_folder):
		os.makedirs(dst_folder)
	for i in range(last_frame+1):
		src=input_folder+"\\"+str(i)+"\\"+target_file
		dst=dst_folder+"\\"+str(i)
		if not os.path.exists(dst):
			os.makedirs(dst)
		copyfile(src, dst+"\\"+target_file)