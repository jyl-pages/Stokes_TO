# usage: python img2gif.py folder fps
# dependencies: PIL package (pip install Pillow)

from PIL import Image, ImageDraw
import sys, glob, time
if __name__ == "__main__":
	folder = sys.argv[1]; fps = int(sys.argv[2]);
	ext = ['png', 'jpg', 'jpeg']
	files = []; [files.extend(glob.glob(folder + '\\*.' + e)) for e in ext]; files.sort();
	gif_file = folder + '\\' + str(int(time.time()))+'.gif'
	images = [Image.open(i) for i in files]
	images[0].save(gif_file, save_all=True, append_images=images[1:], optimize=False, duration=1000/fps)
	print ("Created gif in " + gif_file)