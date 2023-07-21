import imageio
import glob
import os
import sys

def gifs_from_png(list_of_filenames, save_path, duration=0.2):
    images = []

    for filename in list_of_filenames:
        if os.path.exists(filename):
            images.append(imageio.imread(filename))

    imageio.mimsave(save_path, images, duration=duration)

    for f in list_of_filenames:
        if os.path.exists(f):
            os.remove(f)