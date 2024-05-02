import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import argparse
import csv
import subprocess
import platform

from setup import setup
from animate import animate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-save', action='store_true', help='Save animation to file')
    parser.add_argument('test_case', help='Test case to run', choices=['test','uniform','uniform-vertical','vortex','couette','poiseuille','breakup','terminal','custom'])
    args = parser.parse_args()
    test_case = args.test_case
    save_animation = args.save
    
    setup(test_case)
    
    if platform.system() == 'Windows':
        subprocess.run('particles.exe')
    elif platform.system() == 'Darwin':
        subprocess.run('./particles_xcode')
    
    animate('', save_animation)
    return


if __name__ == '__main__':
    main()
