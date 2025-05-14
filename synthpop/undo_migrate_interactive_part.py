"""
This script migrate the model modules constants and data dir to a different location
This might be helpful if synthpop will be installed  at the default package location with pip
e.g. by using pip install.
also the data volume for Isochrones can reach several GB.

It uses symbolic links to connect the two directories.
"""
__all__ = ["undo_migrate"]
__author__ = ["M.J. Huston", "J. Klüter"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__data__ = "2024-10-23"

import sys
import os
import shutil
import readline
from tkinter import Tk, messagebox, TclError
from tkinter.filedialog import askdirectory
from tkinter.simpledialog import Dialog


def copy_dir(src_dir, target_dir, name):
    """copy a directory and create a symlink to the new location"""
    # estimate paths
    src = os.path.join(src_dir, name)
    symlink = src
    des = os.path.join(target_dir, name)

    if des != src:
        # delete link
        os.unlink(des)
        # move directory
        shutil.copytree(src, des, dirs_exist_ok=True)
        shutil.rmtree(src)



def copy_file(src_dir, target_dir, name):
    """copy a file and create a symlink to the new location"""
    # estimate paths
    src = os.path.join(src_dir, name)
    symlink = src
    des = os.path.join(target_dir, name)

    # create a new link
    if des != src:
        # delete link
        os.unlink(des)
        # copy file to new location
        shutil.copy(src, des)
        os.remove(src)


def get_dirname_from_command_line():
    """get the directory from the command line using an autocomplete"""

    def complete(text, state):
        # Search for all files and directories starting with the entered text
        search_dirname, filename = os.path.split(text)
        if search_dirname == '':
            search_dirname2 = '.'
        else:
            search_dirname2 = search_dirname
        options = [os.path.join(search_dirname, f + '/') for f in os.listdir(search_dirname2)
                   if f.startswith(filename)
                   and os.path.isdir(os.path.join(search_dirname2, f))]

        # If there are no options, exit the autocomplete function
        if not options:
            return None

        # If the user presses TAB for the first time, show all options
        if state == 0:
            complete.matches = options

        # Return the next item in the list of options
        if state < len(complete.matches):
            return complete.matches[state]
        else:
            return None

    print("Please specify the directory that holds "
          "the models, modules and configurations, etc.")
    _delims = readline.get_completer_delims()
    readline.parse_and_bind("tab: complete")
    readline.set_completer_delims(' \t\n;')
    readline.set_completer(complete)

    dirname = input('Directory: ')

    readline.set_completer_delims(_delims)
    readline.set_completer(None)
    readline.parse_and_bind("tab: self-insert")

    return dirname


def get_dirname_from_gui():
    """ get filename from a TKinter gui """
    Tk().withdraw()
    result = messagebox.showinfo(
        "Select Directory",
        "Please specify the directory that holds\n"
        "the models, modules and configurations, etc.",
        type=messagebox.OKCANCEL)
    if result == "cancel":
        return
    else:
        dirname = askdirectory(
            initialdir=os.path.expanduser('~'),
            title="test")

    return dirname


def undo_migrate(dirname=''):
    if dirname == '':
        try:
            dirname = get_dirname_from_gui()
        except TclError:
            dirname = get_dirname_from_command_line()

    if dirname == '':
        return

    print(f"Undoing SynthPop migration to {dirname}")
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # copy modules models and config & constants to new location
    synthpop_code_dir = os.path.dirname(__file__)
    copy_dir( dirname,synthpop_code_dir, "config_files")
    copy_dir(dirname, synthpop_code_dir, "modules")
    copy_dir(dirname, synthpop_code_dir, "models")
    copy_file(dirname, synthpop_code_dir, "constants.py")
    os.unlink(synthpop_code_dir+'/outputfiles')
    print("Synthpop_Directory migration has been undone. You can now safely update SynthPop via pip.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        directory = os.path.abspath(sys.argv[1])
    else:
        directory = ''
    undo_migrate(directory)
