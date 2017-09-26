# This is a small script for animating png images. 
# It creates rather large gifs, so you might use other software to make it smaller
# for example, ezgif (a website) 
# You will need to install imageio for this to work
# conda install imageio

import os
import imageio

png_dir = "./examples/spatially_varying_K/"
images = []
for subdir, dirs, files in os.walk(png_dir):
    for file in files:
        file_path = os.path.join(subdir, file)
        if file_path.endswith(".png"):
            images.append(imageio.imread(file_path))
imageio.mimsave('./test.gif', images, subrectangles=True,)
