# MAIN FILE
from tkinter import *
from functools import partial
import numpy as np
import math
import sys
import db

env_dev_flag = False

# BENCENO TEMP VARS
A = 3.98523
B = 1184.24
C = 217.572
Tc = 289.01
Pc = 48.990637

class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

def fenske():
	return 0

def start(f, p, tf, columnaP, xlk, xhk):
	if(not env_dev_flag):
		f = float(f.get())
		p = float(p.get())
		tf = float(tf.get())
		columnaP = float(columnaP.get())
		xlk = float(xlk.get())
		xhk = float(xhk.get())

	fenske()

	# Flujos de alimentacion de componentes ligero y pesado
	xlf = 0.5
	xhf = 0.8
	flf = xlf * f
	fhf = xhf * f

	# Flujos de componentes en el domo y fondo
	fld = flf * xlk
	fhd = fhf * xhk
	flw = flf - fld
	fhw = fhf - fhd

	temperaturaBurbuja(p, zf, tf)
	return


def temperaturaBurbuja(p, zi, tf):
	# Fijar P y composición
	xi = zi
	# Estimar Td
	td = tf
	# Calcular Ki como primera aprox
	pi_sat = math.exp(A - (B / (tf + C)))
	ki = pi_sat / p
	# Calcular yi_supuesta
	yi_supuesta = ki * xi
	
	while(True):
		kib = calculateKib(td, Tc, Pc, xi, yi_supuesta, p)
		break

def calculateKib(t, tci, pci, xi, yi_supuesta, p):
	print(t)
	Ai = (0.4278/(pci*((t/tci)**2.5)))**(1/2)
	Bi = 0.0867/(pci*(t/tci))
	print("Ai, Bi")
	print(Ai, Bi)
	Bx = 0.0
	By = 0.0
	Ax = 0.0
	Ay = 0.0
	#for i in range(0, int(C)): # C es de componentes elegidos
	for i in range(0, 4): # C es de componentes elegidos
		Bx += Bi * xi
		By += Bi * yi_supuesta
		Ax += Ai * xi
		Ay += Ai * yi_supuesta
	
	this_A = math.sqrt(Ax * Ay)
	this_B = math.sqrt(Bx * By)
	coef_c = this_B * p * ( (this_A**2/this_B) - (this_B*p) - 1)
	coef_d = -1 * (this_A**2/this_B) * (this_B*p)**2
	print(coef_c, coef_d)
	print("Roots of the polynomial:")
	roots = np.roots([1, -1, coef_c, coef_d])
	print(roots)
	for i in range(0, len(roots)):
		print("Root:", float(roots[i]))
	return 0

def print_selection(key):
	var_obj = var.get(key)
	print(var_obj.get())

def main():
    print('Bienvenido: Simulador')
    db.open()
    elements = db.getElements()
    if(env_dev_flag):
    	db.open()
    	db.getElements()
    	f = 100
    	zf = 0.8
    	p = 102
    	tf = 273.15
    	columnaP = 1
    	xlk = 5.2
    	xhk = 2.5
    	start(f, zf, p, tf, columnaP, tipoCondensador, xlk, xhk)
    	sys.exit()
	
    root = Tk()
    ##app = Window(root)
    root.wm_title("Simulador") # set window title

    # Define form
    FLabel = Label(root, text="F").grid(row=0, column=0)
    f = StringVar()
    usernameEntry = Entry(root, textvariable=f).grid(row=0, column=1)

    zfLabel = Label(root,text="Zf").grid(row=1, column=0)
    zfLabel = Label(root,text="XLF").grid(row=1, column=3)
    zfLabel = Label(root,text="XHF").grid(row=1, column=4)

    selectCheck = dict()
    xlfCheck = dict()
    xhfCheck = dict()
    zfDicc = dict()
    zfInputList = list()
    rowCounter = 2
    elementDicc = dict()
    for i in range(len(elements)):
        rowCounter += i
        element = elements[i]
        selectCheck[element] = IntVar()
        xlfCheck[element] = IntVar()
        xhfCheck[element] = IntVar()
        zfDicc[element] = StringVar()
        Checkbutton(root, variable=selectCheck[element]).grid(row=rowCounter, column=0)
        Label(root, text=str(element)).grid(row=rowCounter, column=1)
        zfEntry = Entry(root, textvariable=zfDicc[element]).grid(row=rowCounter, column=2) 
        Checkbutton(root, variable=xlfCheck[element]).grid(row=rowCounter, column=3)
        Checkbutton(root, variable=xhfCheck[element]).grid(row=rowCounter, column=4)
        zfInputList.append(zfEntry)
    
    rowCounter += 1
    pLabel = Label(root,text="P").grid(row=rowCounter, column=0)  
    p = StringVar()
    pEntry = Entry(root, textvariable=p).grid(row=rowCounter, column=1) 

    rowCounter += 1
    tfLabel = Label(root,text="Tf").grid(row=rowCounter, column=0)  
    tf = StringVar()
    tfEntry = Entry(root, textvariable=tf).grid(row=rowCounter, column=1)

    rowCounter += 1
    columnaPLabel = Label(root,text="Columna P").grid(row=rowCounter, column=0)  
    columnaP = StringVar()
    columnaPEntry = Entry(root, textvariable=columnaP).grid(row=rowCounter, column=1) 

    rowCounter += 1
    xlkLabel = Label(root,text="XLK").grid(row=rowCounter, column=0)  
    xlk = StringVar()
    xlkEntry = Entry(root, textvariable=xlk).grid(row=rowCounter, column=1) 

    rowCounter += 1
    xhkLabel = Label(root,text="XHK").grid(row=rowCounter, column=0)  
    xhk = StringVar()
    xhkEntry = Entry(root, textvariable=xhk).grid(row=rowCounter, column=1)

    rowCounter += 1
    ropLabel = Label(root,text="R OP").grid(row=rowCounter, column=0)  
    r_op = StringVar()
    xhkEntry = Entry(root, textvariable=r_op).grid(row=rowCounter, column=1) 

    startsSimulation = partial(start, f, p, tf, columnaP, xlk, xhk)

    rowCounter += 1
    startButton = Button(root, text="Simular", 
    	highlightbackground='#000000', 
    	command=startsSimulation
	).grid(row=rowCounter, column=1)  
    #show window
    print("Show window")
    root.mainloop()
    print('Programa terminado')
 
if __name__ == "__main__":
    main()