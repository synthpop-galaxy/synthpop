"""
This script compares a default config file specified by the user
with the _default.synthpop_conf file and highlights missing items.
"""
import sys
import os
import json
import argparse

# get the _default_file
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT = os.path.join(THIS_DIR, "_default.synthpop_conf")

def compare_json_files(userfile_path):
    with open(_DEFAULT, 'r') as f:
        default_file = json.load(f)
    with open(userfile_path, 'r') as f:
        userfile = json.load(f)
    is_ok = True
    for section in default_file.keys():
        if section.startswith('#'): continue
        default_keys = {item for item in default_file[section].keys() if not item.startswith('#')}
        # check if section is in userfile
        if section not in userfile.keys():
            is_ok = False
            print(f"'{userfile_path}' should contain a '{section}' section!")
            print(f"this should included the following attributes:")
            print(default_keys)
            continue
        user_key = {item for item in userfile[section].keys() if not item.startswith('#')}
        missing_attributes = default_keys - user_key
        if missing_attributes:
            is_ok=False
            print(f"'{section}' section must contain the following attributes:")
            print(missing_attributes)
        extra_attributes = user_key - default_keys
        if extra_attributes:
            print(f"'{section}' section contain the following extra attributes:")
            print(extra_attributes)
            print(f"these might lead to an unexpected behavior.")

    extra_sections = set(userfile)-    set(default_file.keys())
    if extra_sections:
        print(f"'{userfile_path}' contains the following extra sections:")
        print(extra_sections)
        print(f"these might lead to an unexpected behavior.")

    if is_ok:
        print(f"\n'{userfile_path}' can be used as default_config file")
    else:
        print('\nWARNING!')
        print(f"{userfile_path}' should NOT be used as default_config file!")

    return is_ok


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Compares a config file with _default.synthpop_conf file and determines '
                    'if it can be used as default config file ')

    # add the argument for userfile_path
    parser.add_argument('userfile_path', type=str, help='the path to the file')

    # parse the arguments
    args = parser.parse_args()

    compare_json_files(args.userfile_path)
