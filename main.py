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

    print('-> Calcular temp de burbuja')
    it_count = 0
    result_t_burbuja = temperaturaBurbuja(p, zfDicc, tf)
    while not result_t_burbuja['status']:
        it_count+=1
        print('Iteracion:', it_count)
        result_t_burbuja = temperaturaBurbuja(p, zfDicc, result_t_burbuja['Td'])

    print('     Regrese pero no hay nadie:', Td_obtenida)
    input()
    return

def temperaturaBurbuja(p, zfDicc, tf):
    Td_obtenida = 0.0
    yis_calculadas = list()
    cicle_counter = 1
    result = dict()
    for element in zfDicc:
        # Fijar P y composición
        xi = float(zfDicc[element]['value'].get())
        # Estimar Td
        td = tf
        # Obtener valores de BD para el elemento
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        # Calcular Ki como primera aprox
        Ki = calculate_Ki(p, td, db_values)
        # Calcular yi_supuesta
        yi_supuesta = calculate_yi_supuesta(Ki, xi)
        # Obtner Kib
        Kib = calculate_Kib(db_values, zfDicc, td, p, yi_supuesta)
        while True:
            # Calcular yi_calculada
            yi_calculada=Kib*xi
            #print('yi calculada = Kib·xi')
            #print('yi calculada = ', Kib,'*',xi, '=', yi_calculada)
            normalizado = abs(yi_supuesta - yi_calculada)
            print('abs(yi_supuesta - yi_calculada) => ', yi_supuesta, '-', yi_calculada, '=', normalizado)
            if(normalizado  < 0.001):
                #print('<<-- Yi calculada para', element, 'en ciclo', cicle_counter)
                #input()
                yis_calculadas.append(yi_calculada)
                #print('->Yis calculadas totales:', len(yis_calculadas), 'Cantidad de componentes:', len(zfDicc))
                if len(yis_calculadas) == len(zfDicc):
                    f = 1 - sum(yis_calculadas)
                    if f >= -0.01 and f <= 0.01:
                        print('<------------------------------------------------>')
                        Td_obtenida = td
                        result['status'] = True
                        result['Td'] = Td_obtenida
                        result['yis_calculadas'] = yis_calculadas
                        return result
                    else:
                        delta_td = 0
                        if f < 0:
                            print('delta_td -0.01')
                            delta_td = -0.01
                        else:
                            print('delta_td +0.01')
                            delta_td = 0.01
                        td = td + delta_td
                        print('Td=', td)
                        result['status'] = False
                        result['Td'] = td
                        return result
                else:
                    break
            else:
                yi_supuesta = yi_calculada
                Kib = calculate_Kib(db_values, zfDicc, td, p, yi_supuesta)
                print(Kib)
                cicle_counter+=1
                print('Ciclo', cicle_counter, 'falló')
                if(cicle_counter == 10): sys.exit()
    print('Proceso recursivo ha fallado.')
    sys.exit()

def calculate_Ki(p, td, db_values):
    A = db_values['A']
    B = db_values['B']
    C = db_values['C']
    pi_sat = math.exp(A - (B / (td + C)))
    Ki = pi_sat / p
    return Ki

def calculate_yi_supuesta(ki, xi):
    yi_supuesta = ki * xi
    #print('yi_supuesta = Ki·xi')
    #print('yi_supuesta = ', ki,'*',xi, '=', yi_supuesta)
    return yi_supuesta

def calculate_Kib(db_values, zfDicc, tf, p, yi_supuesta):
    Tc = db_values['Tc']
    Pc = db_values['Pc']
    Ai = calculate_Ai(Pc, tf, Tc)
    Ay = calculate_Ay(zfDicc, yi_supuesta, tf)
    Ax = calculate_Ax(zfDicc, tf)
    Bi = calculate_Bi(Pc, tf, Tc)
    By = calculate_By(zfDicc, yi_supuesta, tf)
    Bx = calculate_Bx(zfDicc, tf)

    A = math.sqrt(Ax * Ay)
    B = math.sqrt(Bx * By)

    Z = getRoots(A, B, p)
    Zl = Z['Zl']
    Zv = Z['Zv']
    FIv = coeficiente_de_fugacidad(Zv, p, A, B, Ai, Bi)
    FIl = coeficiente_de_fugacidad(Zl, p, A, B, Ai, Bi)
    print('FIv:', FIv)
    print('FIl:', FIl)
    return FIl/FIv

def calculate_Ai(Pc, tf, Tc):
    return (0.4278/(Pc*((tf/Tc)**2.5)))**(0.5)

def calculate_Bi(Pc, tf, Tc):
    return 0.0867/(Pc*(tf/Tc))

def calculate_Bx(zfDicc, tf):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Bi = calculate_Bi(Pc, tf, Tc)
        sumatoria += Bi * float(zfDicc[element]['value'].get())
    return sumatoria

def calculate_Ax(zfDicc, tf):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Ai = calculate_Ai(Pc, tf, Tc)
        sumatoria += Ai * float(zfDicc[element]['value'].get())
    return sumatoria

def calculate_Ay(zfDicc, yi_supuesta, tf):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Ai = calculate_Ai(Pc, tf, Tc)
        sumatoria += Ai * yi_supuesta
    return sumatoria

def calculate_By(zfDicc, yi_supuesta, tf):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Bi = calculate_Bi(Pc, tf, Tc)
        sumatoria += Bi * yi_supuesta
    return sumatoria

def getRoots(A, B, p):
    coef_c = (B * p) * ( (A**2/B) - (B*p) - 1)
    coef_d = -1 * ((A**2/B) * ((B*p)**2))
    roots = np.roots([1, -1, coef_c, coef_d])
    # for i in range(0, len(roots)):
    #     roots[i] = float(roots[i])
    roots.sort()
    print('Raices:', roots)
    Z = dict()
    Z['Zl'] = roots[2]
    Z['Zv'] = roots[1]
    print('Zv:',Z['Zv'])
    print('Zl:',Z['Zl'])
    return Z

def coeficiente_de_fugacidad(Z, p, A, B, Ai, Bi):
    first_eq = (Z - 1) * (Bi / B)
    second_eq = math.log(Z - (B * p))
    subeq3_1 = (A**2 / B) * ((2*Ai / A) - (Bi / B))
    subeq3_2 = math.log(1 + ((B * p) / Z))
    third_eq = subeq3_1 * subeq3_2
    final_eq = first_eq - second_eq - third_eq
    return math.exp(final_eq)

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
        f.set(100)
        for element in elementDicc:
            elementDicc[element]['is_select'].set(1)
            if element == 'n-Hexano':
                print('----->   Set Etano to XLF',)
                elementDicc[element]['is_xlf'].set(1)
            if element == 'o-Xileno':
                print('----->   Set Cumeno to XHF',)
                elementDicc[element]['is_xhf'].set(1)

        elementDicc['Benceno']['value'].set(0.05)
        elementDicc['Acetona']['value'].set(0.01)
        elementDicc['n-Butano']['value'].set(0.15)
        elementDicc['Etano']['value'].set(0.05)
        elementDicc['Cumeno']['value'].set(0.1)
        elementDicc['n-Heptano']['value'].set(0.1)
        elementDicc['n-Hexano']['value'].set(0.1)
        elementDicc['o-Xileno']['value'].set(0.1)
        elementDicc['Acetato de Etilo']['value'].set(0.1)
        p.set(1)
        tf.set(49.8)
        columnaP.set(1)
        xlk.set(0.97)
        xhk.set(0.03)
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