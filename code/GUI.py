import tkinter
from tkinter import *
from tkinter import filedialog
import crystal_visual
import pandas as pd
import numpy as np
import os
from tkinter import font as tkFont

window = tkinter.Tk()
window.title("Crystal Visual")
window.geometry("600x800")

file_path = tkinter.Entry(window)
# file_path.place(x = 20, y = 20)

sys20 = tkFont.Font(family='system', size=15, weight=tkFont.BOLD)
sys15_nor = tkFont.Font(family='system', size=15, weight=tkFont.NORMAL)
sys36 = tkFont.Font(family='system', size=36, weight=tkFont.BOLD)
sys25 = tkFont.Font(family='system', size=25, weight=tkFont.BOLD)

tkinter.Label(window, text = 'Number of expanding', font = sys20).place(x=220, y=100)
option = [num for num in range(3,20) if num%2==1]
variable_num_expanding = StringVar(window)
variable_num_expanding.set(option[0])
num_expanding = OptionMenu(window, variable_num_expanding, *option)
num_expanding.place(x= 400, y = 100)


tkinter.Label(window, text = 'Number of saved distance', font = sys20).place(x=220, y=160)
option = range(1,201)
variable_num_dist = StringVar(window)
variable_num_dist.set(option[0])
num_dist = OptionMenu(window, variable_num_dist, *option)
num_dist.place(x= 430, y = 160)

tkinter.Label(window, text = 'Type of atom', font = sys20).place(x=220, y=220)
type_atom = tkinter.Entry(window)
type_atom.place(x = 220, y=260)

tkinter.Label(window, text = 'Files', font = sys15_nor).place(x=10, y=80)
fileExploer = tkinter.Listbox(window, height = 30, width = 20)
fileExploer.place(x = 10, y = 100)


tkinter.Button(window, text = 'Select folder', font = sys36, command = lambda:open_folder()).place(x = 10, y = 25)
tkinter.Button(window, text = 'Crystal Visulization', font = sys25, command = lambda:crystal_visulization()).place(x = 220, y = 400)
tkinter.Button(window, text = 'Cluster Crystals', font = sys25, command = lambda:classify_crystal()).place(x = 220, y =500)



def crystal_visulization():
    files_name = os.path.join(file_path.get(),fileExploer.get(fileExploer.curselection()))
    num_dist = int(variable_num_dist.get())
    num_expand=int(variable_num_expanding.get())
    atom_type = type_atom.get()
    num_file = 1
    if len(atom_type) == 0:
        atom_type = None
        crystal_visual.final_step(files = files_name, 
        atom_type = atom_type,num_file = num_file, num_dist=num_dist, num_expand=num_expand)
    else:
        crystal_visual.final_step(files = files_name, 
        atom_type = [x.strip() for x in atom_type.split(',')], num_file = num_file, num_dist=num_dist, num_expand=num_expand)

def classify_crystal():
    path = file_path.get()
    num_dist = int(variable_num_dist.get())
    num_expand=int(variable_num_expanding.get())
    atom_type = type_atom.get()
    num_file = fileExploer.size()
    if len(atom_type) == 0:
        atom_type = None
        crystal_visual.final_step(files = path, 
        atom_type = atom_type, num_dist=num_dist, num_expand=num_expand, num_file = num_file)
    else:
        crystal_visual.final_step(files = path, 
        atom_type = [x.strip() for x in atom_type.split(',')], num_dist=num_dist, num_expand=num_expand, num_file = num_file)



def open_folder():
    file_select =  filedialog.askdirectory()
    file_path.insert(tkinter.END, file_select)
    show_dir(file_path.get())
 

def show_dir(path):
        i = 0
        for r, d, f in os.walk(path):
            for file in f:
                if file.endswith('.cif'):
                    i = i+1
                    fileExploer.insert(i, file) 

window.mainloop()
