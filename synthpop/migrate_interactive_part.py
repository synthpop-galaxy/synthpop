"""
This script migrate the model modules constants and data dir to a different location
This might be helpful if synthpop will be installed  at the default package location with pip
eg by using pip install.
also the data volumne for Isochrones can reach several GB.

It uses symbolic links to connect the two directories.
"""


import sys
import os
import shutil 
from tkinter import Tk, messagebox
from tkinter.filedialog import askdirectory
from tkinter.simpledialog import Dialog

def copydir(src_dir, target_dir, name):
    """copy a directory and create a symlik to the new location"""
    # estimate paths
    src = os.path.join(src_dir, '_' + name)
    symlink = os.path.join(src_dir, name)
    des = os.path.join(target_dir, name)

    # deleat previous link
    if os.path.islink(symlink):
        os.unlink(symlink)

    if des != symlink:
        # copy directory to new location
        shutil.copytree(src, des, dirs_exist_ok = True)
        # create a new link
        os.symlink(des, symlink, target_is_directory = True)
    else:
        shutil.move(src, des)

def copyfile(src_dir, target_dir, name):
    """copy a file and create a symlik to the new location"""
    # estimate paths
    src = os.path.join(src_dir, '_' + name)
    symlink = os.path.join(src_dir, name)
    des = os.path.join(target_dir, name)

    # deleat previous link
    if os.path.islink(symlink):
        os.unlink(symlink)

    if des != symlink:
        # copy file to new location
        shutil.copy(src, des)
        # create a new link

        os.symlink(des, symlink, target_is_directory = True)
    else:
        shutil.move(src, des)
def migrate(dirname = None):
    if dirname is None:
        # get filename from gui
        Tk().withdraw()
        result = messagebox.showinfo("Select Directory",
                "Please specify a directory for easy axcess "
                "to the models, modules and configurations, etc.",
                type=messagebox.OKCANCEL)
        if result == "cancel":
            return
        else:
            dirname = askdirectory(
            initialdir=os.path.expanduser('~'),
            title="test")
    if dirname == '':
        return
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # copy modules models and config & constants to new location
    synthpop_code_dir = os.path.dirname(__file__)
    copydir(synthpop_code_dir, dirname, "config_files")
    copydir(synthpop_code_dir, dirname, "modules")
    copydir(synthpop_code_dir, dirname, "models")
    copyfile(synthpop_code_dir, dirname, "constants.py")

if __name__=="__main__":
    if len(sys.argv) > 1:
        dirname = sys.argv[1]
    else: 
        dirname = None
    migrate(dirname)

