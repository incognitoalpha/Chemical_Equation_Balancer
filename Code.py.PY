import re
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk
from sympy import Matrix, lcm
import matplotlib.pyplot as plt
import csv


try:
    element_data = pd.read_csv('PubChemElements_all.csv')
except FileNotFoundError:
    print("Error: PubChemElements_all.csv file not found.")
    exit()
except pd.errors.EmptyDataError:
    print("Error: PubChemElements_all.csv file is empty or invalid.")
    exit()

symbolList = []
symbolMatrix = []

def set_dark_theme():
    style = ttk.Style()
    style.configure("TButton", foreground="black", background="black")
    style.configure("TEntry", fieldbackground="#222222", foreground="BLACK")
    style.configure("TLabel", foreground="white", background="black")

def addToMatrix(symbol, index, count, side):
    if(index == len(symbolMatrix)):
        symbolMatrix.append([])
        for x in symbolList:
            symbolMatrix[index].append(0)
    if(symbol not in symbolList):
        symbolList.append(symbol)
        for i in range(len(symbolMatrix)):
            symbolMatrix[i].append(0)
    column = symbolList.index(symbol)
    symbolMatrix[index][column] += count * side

def findSymbols(segment, index, multiplier, side):
    elementsAndNumbers = re.split('([A-Z][a-z]?)', segment)
    i = 0
    while(i < len(elementsAndNumbers)-1):
        i += 1
        if(len(elementsAndNumbers[i]) > 0):
            if(elementsAndNumbers[i+1].isdigit()):
                count = int(elementsAndNumbers[i+1]) * multiplier
                addToMatrix(elementsAndNumbers[i], index, count, side)
                i += 1
            else:
                addToMatrix(elementsAndNumbers[i], index, multiplier, side)

def compoundDecipher(compound, index, side):
    segments = re.split(r'(\([A-Za-z0-9]*\)[0-9]*)', re.sub(r'^\d+', '', compound))
    for segment in segments:
        if segment.startswith("("):
            segment = re.split(r'\)([0-9]*)', segment)
            multiplier = int(segment[1])
            segment = segment[0][1:]
        else:
            multiplier = 1
        findSymbols(segment, index, multiplier, side)

class ChemicalApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemical Equation and Element Properties")
        self.root.geometry("1000x600")
        self.root.resizable(False, False)

        self.canvas = tk.Canvas(root, width=1000, height=600)
        self.canvas.configure(background='black')
        self.canvas.pack()

        try:
            background_image = tk.PhotoImage(file="istockphoto-1318101781-170667a.png")
            self.canvas.create_image(0, 0, anchor=tk.NW, image=background_image)
        except tk.TclError:
            print("Error: Image file not found or unsupported format.")
            exit()

        reactants_label = ttk.Label(root, text="Enter Reactants (In skeletal format):", font=("Arial", 12))
        reactants_label.place(x=10, y=10)

        self.reactants_entry_var = tk.StringVar()
        self.reactants_entry = ttk.Entry(root, textvariable=self.reactants_entry_var, font=("Arial", 12))
        self.reactants_entry.place(x=270, y=10)

        products_label = ttk.Label(root, text="Enter Products (In skeletal format):", font=("Arial", 12))
        products_label.place(x=10, y=50)

        self.products_entry_var = tk.StringVar()
        self.products_entry = ttk.Entry(root, textvariable=self.products_entry_var, font=("Arial", 12))
        self.products_entry.place(x=270, y=50)

        submit_button = ttk.Button(root, text="Submit", command=self.process_equation, style="TButton")
        submit_button.place(x=220, y=90)

        self.equation_label = ttk.Label(root, text="", font=("Arial", 14, "bold"))
        self.equation_label.place(x=10, y=130)

        output_properties = ttk.Label(root, text="Properties of Unique Elements:", font=("Arial", 14, "underline"))
        output_properties.place(x=10, y=170)

        scrollbar_x = ttk.Scrollbar(root, orient="horizontal")

        self.elements_table = ttk.Treeview(root, xscrollcommand=scrollbar_x.set)
        self.elements_table['columns'] = [col for col in element_data.columns if col != 'CPKHexColor']

        for column in self.elements_table['columns']:
            self.elements_table.heading(column, text=column, anchor='center')
            self.elements_table.column(column, anchor='center', width=100)

        scrollbar_x.config(command=self.elements_table.xview)
        scrollbar_x.place(x=10, y=210, width=980)

        self.elements_table.place(x=10, y=240, width=980)

    def process_equation(self):
        symbolList.clear()
        global symbolMatrix
        symbolMatrix = []

        for item in self.elements_table.get_children():
            self.elements_table.delete(item)
        self.equation_label.config(text="")

        reactants = [re.sub(r'^\d+', '', string) for string in self.reactants_entry_var.get().replace(' ', '').split("+")]
        products = [re.sub(r'^\d+', '', string) for string in self.products_entry_var.get().replace(' ', '').split("+")]

        for i in range(len(reactants)):
            compoundDecipher(reactants[i], i, 1)

        for i in range(len(products)):
            compoundDecipher(products[i], i + len(reactants), -1)

        symbolMatrix = Matrix(symbolMatrix)
        symbolMatrix = symbolMatrix.transpose()
        solution = symbolMatrix.nullspace()[0]
        multiple = lcm([val.q for val in solution])
        solution = multiple * solution
        coEffi = solution.tolist()

        equation_text = ""
        for i in range(len(reactants)):
            equation_text += str(coEffi[i][0]) + reactants[i]
            if i < len(reactants) - 1:
                equation_text += " + "
        equation_text += " -> "
        for i in range(len(products)):
            equation_text += str(coEffi[i + len(reactants)][0]) + products[i]
            if i < len(products) - 1:
                equation_text += " + "
    
        self.equation_label.config(text=equation_text)

        unique_symbols = set(symbolList)

        for symbol in unique_symbols:
            symbol_info = element_data[element_data['Symbol'] == symbol]
            if not symbol_info.empty:
                properties = symbol_info.iloc[0].to_dict()

                if 'CPKHexColor' in properties:
                    properties.pop('CPKHexColor')

                self.elements_table.insert('', 'end', text=symbol,
                                           values=[properties.get(col, '') for col in self.elements_table['columns']])

        self.plot_graphs()

    def plot_graphs(self):       
        unique_symbols_data = []
        symbol_names = []
        
        
        for symbol in self.elements_table.get_children():
            #print("symbol",symbol)
            values = self.elements_table.item(symbol, 'values')
            # print("values",values)
            unique_symbols_data.append({self.elements_table['columns'][i]: values[i] for i in range(len(values))})
            # print( "unique_symbol_data",unique_symbols_data)
            symbol_names.append(values[1])
            # print("symbol_names",symbol_names)
            
            
        properties_to_plot = ['AtomicNumber', 'AtomicMass', 'Electronegativity', 'AtomicRadius',
                          'IonizationEnergy', 'ElectronAffinity', 'MeltingPoint', 'BoilingPoint', 'Density']
        property_labels = ['Atomic Number', 'Atomic Mass', 'Electronegativity', 'Atomic Radius',
                       'Ionization Energy', 'Electron Affinity', 'Melting Point', 'Boiling Point', 'Density']
        
        
        
        fig, axs = plt.subplots(3, 3, figsize=(15, 12))
        
        
        
        for i, prop in enumerate(properties_to_plot):
            row = i // 3
            col = i % 3
            ax = axs[row, col]
            ax.set_title(property_labels[i] + ' Variation')


           

            sorted_data = sorted(zip(symbol_names, [data.get(prop, 0) for data in unique_symbols_data]), key=lambda x: x[1])
            
            sorted_names, sorted_values = zip(*sorted_data)
            


            
               
            ax.scatter(sorted_names, sorted_values, label=property_labels[i])
            
            
            
            ax.set_xlabel('')
            ax.set_ylabel(property_labels[i])
            ax.tick_params(axis='x', rotation=45)
            ax.legend()
            
            
        plt.tight_layout()
        plt.show()       

root = tk.Tk()
set_dark_theme()
app = ChemicalApp(root)
root.mainloop()
