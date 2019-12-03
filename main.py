# MAIN FILE
from tkinter import *
from functools import partial
import numpy as np
import math
import sys
import db

env_dev_flag = True

class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

def purgeDicc(zfDicc):
    r = dict(zfDicc)
    for item in zfDicc:
        if(not zfDicc[item]['is_select'].get() == 1):
            del r[item]
    return r

def start(f, zfDicc, p, tf, columnaP, xlk, xhk):

    zfDicc = purgeDicc(zfDicc)
    f = float(f.get())
    p = float(p.get())
    tf = float(tf.get())
    columnaP = float(columnaP.get())
    xlk = float(xlk.get())
    xhk = float(xhk.get())

    # Flujos de alimentacion de componentes ligero y pesado
    xlf = 0.0
    xhf = 0.0
    for item in zfDicc:
        if zfDicc[item]['is_xlf'].get() == 1:
            xlf = float(zfDicc[item]['value'].get())
        if zfDicc[item]['is_xhf'].get() == 1:
            xhf = float(zfDicc[item]['value'].get())
    
    flf = xlf * f
    fhf = xhf * f

	# Flujos de componentes en el domo y fondo
    fld = flf * xlk
    fhd = fhf * xhk
    flw = flf - fld
    fhw = fhf - fhd

    temperaturaBurbuja(p, zfDicc, tf)
    return


def temperaturaBurbuja(p, zfDicc, tf):
    Bx = 0.0
    By = 0.0
    Ax = 0.0
    Ay = 0.0
    Ai = 0.0
    Bi = 0.0
    ki = 0.0
    for element in zfDicc:
        print('Working with:', element, '=', zfDicc[element]['value'].get())
        # Fijar P y composición
        xi = float(zfDicc[element]['value'].get())
        # Estimar Td
        td = tf
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        A = db_values['A']
        B = db_values['B']
        C = db_values['C']
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        # Calcular Ki como primera aprox
        pi_sat = math.exp(A - (B / (tf + C)))
        ki = pi_sat / p
        # Calcular yi_supuesta
        yi_supuesta = ki * xi
        print('yi_supuesta = Ki·xi')
        print('yi_supuesta = ', ki,'*',xi, '=', yi_supuesta)

        Ai = (0.4278/(Pc*((tf/Tc)**2.5)))**(1/2)
        Bi = 0.0867/(Pc*(tf/Tc))

        Bx += Bi * xi
        By += Bi * yi_supuesta
        Ax += Ai * xi
        Ay += Ai * yi_supuesta


    this_A = math.sqrt(Ax * Ay)
    this_B = math.sqrt(Bx * By)
    coef_c = (this_B * p) * ( (this_A**2/this_B) - (this_B*p) - 1)
    #coef_d = -1 * ((this_A**2/this_B) * ((this_B*p)**2))
    coef_d = -1 * ((this_A**2/this_B) * ((this_B*p)**2))
    roots = np.roots([1, -1, coef_c, coef_d])
    roots.sort()
    print('Raices:', roots)
    Zl = roots[0]
    Zv = roots[1]
    print('Zv',Zv)
    print('Zl',Zl)
    FIv = coeficiente_de_fugacidad(Zv, p, this_A, this_B, Ai, Bi)
    # FIv = math.exp(((Zv - 1) * (Bi / this_B)) - math.log(Zv - (this_B * p)) - ((this_A**2 / this_B) * ((2*Ai / this_A) - (Bi / this_B)) * (math.log(1 + ((this_B * p) / Zv)))))
    FIl = coeficiente_de_fugacidad(Zl, p, this_A, this_B, Ai, Bi)
    # FIl = math.exp(((Zl - 1) * (Bi / this_B)) - math.log(Zl - (this_B * p)) - ((this_A**2 / this_B) * ((2*Ai / this_A) - (Bi / this_B)) * (math.log(1 + ((this_B * p) / Zl)))))
    print('FIv:', FIv, 'FIl:', FIl)
    Kib = FIl/FIv
    yi_calculada=ki*xi
    print('yi calculada = Ki·xi')
    print('yi calculada = ', ki,'*',xi, '=', yi_calculada)
    input()

def coeficiente_de_fugacidad(Z, p, this_A, this_B, Ai, Bi):    
    return math.exp(((Z - 1) * (Bi / this_B)) - math.log(Z - (this_B * p)) - ((this_A**2 / this_B) * ((2*Ai / this_A) - (Bi / this_B)) * (math.log(1 + ((this_B * p) / Z)))))

def main():
    print('Bienvenido: Simulador')
    db.open()
    elements = db.getElements()
	
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

    rowCounter = 2
    elementDicc = dict()
    for i in range(len(elements)):
        rowCounter += i
        element = elements[i]
        elementDicc[element] = dict()
        elementDicc[element]['is_select'] = IntVar()
        elementDicc[element]['is_xlf'] = IntVar()
        elementDicc[element]['is_xhf'] = IntVar()
        elementDicc[element]['value'] = StringVar()
        elementDicc[element]['db_row'] = i+1

        Checkbutton(root, variable=elementDicc[element]['is_select']).grid(row=rowCounter, column=0)
        Label(root, text=str(element)).grid(row=rowCounter, column=1)
        zfEntry = Entry(root, textvariable=elementDicc[element]['value']).grid(row=rowCounter, column=2) 
        Checkbutton(root, variable=elementDicc[element]['is_xlf']).grid(row=rowCounter, column=3)
        Checkbutton(root, variable=elementDicc[element]['is_xhf']).grid(row=rowCounter, column=4)
    
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

    if(env_dev_flag):
        f.set(0.5)
        elementDicc['Benceno']['is_select'].set(1)
        elementDicc['Acetona']['is_select'].set(1)
        elementDicc['Etano']['is_select'].set(1)
        elementDicc['Benceno']['is_xhf'].set(1)
        elementDicc['Etano']['is_xlf'].set(1)
        elementDicc['Benceno']['value'].set(0.1)
        elementDicc['Acetona']['value'].set(0.3)
        elementDicc['Etano']['value'].set(0.5)
        p.set(0.6)
        tf.set(0.5)
        columnaP.set(0.3)
        xlk.set(0.8)
        xhk.set(0.7)
    startsSimulation = partial(start, f, elementDicc, p, tf, columnaP, xlk, xhk)

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