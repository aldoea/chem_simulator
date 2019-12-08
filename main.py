# MAIN FILE
from tkinter import *
from functools import partial
from burbuja import temperaturaBurbuja
#import sys
import db
from pprint import pprint

env_dev_flag = True
composiciones_fondo = dict()
composiciones_fondo['Benceno'] = 0.03
composiciones_fondo['Acetona'] = 0
composiciones_fondo['n-Butano'] = 0
composiciones_fondo['Etano'] = 0
composiciones_fondo['Cumeno'] = 0.19
composiciones_fondo['n-Heptano'] = 0.05
composiciones_fondo['n-Hexano'] = 0.01
composiciones_fondo['o-Xileno'] = 0.37
composiciones_fondo['Acetato de Etilo'] = 0.02
composiciones_fondo['n-Octano'] = 0.33

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
    xlf_name = None
    xhf = 0.0
    xhf_name = None
    for element in zfDicc:
        if zfDicc[element]['is_xlf'].get() == 1:
            xlf = float(zfDicc[element]['value'].get())
            xlf_name = element
        if zfDicc[element]['is_xhf'].get() == 1:
            xhf = float(zfDicc[element]['value'].get())
            xhf_name = element
    
    flf = xlf * f
    fhf = xhf * f

	# Flujos de componentes en el domo y fondo
    fld = flf * xlk
    fhd = fhf * xhk
    flw = flf - fld
    fhw = fhf - fhd

    print('-> Calcular temp de burbuja')
    it_count = 0
    result_t_burbuja = temperaturaBurbuja(columnaP, zfDicc, tf)
    while not result_t_burbuja['status']:
        it_count+=1
        print('Iteracion:', it_count)
        result_t_burbuja = temperaturaBurbuja(columnaP, zfDicc, result_t_burbuja['Td'])
    print('TD Obtenida:', result_t_burbuja['Td'])
    
    Td_obtenida = result_t_burbuja['Td']
    #yi_obtenida = result_t_burbuja['yis_calculadas']
    kibs_d_obtenida = result_t_burbuja['kib_calculadas']
    
    # Copia de diccionario
    it_count = 0
    dict_componentes = dict(zfDicc)
    for element in dict_componentes:
        xiw = composiciones_fondo[element]
        dict_componentes[element]['value'] = StringVar()
        dict_componentes[element]['value'].set(xiw)

    result_t_burbuja = temperaturaBurbuja(columnaP, dict_componentes, tf)
    while not result_t_burbuja['status']:
        it_count+=1
        print('Iteracion:', it_count)
        result_t_burbuja = temperaturaBurbuja(columnaP, dict_componentes, result_t_burbuja['Td'])
    print('TW Obtenida:', result_t_burbuja['Td'])

    Tw_obtenida = result_t_burbuja['Td']
    kibs_w_obtenida = result_t_burbuja['kib_calculadas']

    alfa_iD = 0.0
    alfa_iW = 0.0
    print('Kibs D:')
    pprint(kibs_d_obtenida)
    print('Kibs W:')
    pprint(kibs_w_obtenida)
    print('Pesado:', xhf_name)
    KidH = kibs_d_obtenida[xhf_name]
    KiwH = kibs_w_obtenida[xhf_name]
    print(KidH, KiwH)
    for element in zfDicc:
        alfa_iD = kibs_d_obtenida[element] / KidH
        alfa_iW = kibs_w_obtenida[element] / KiwH
        zfDicc[element]['alfa_iD'] = alfa_iD
        zfDicc[element]['alfa_iW'] = alfa_iW
        print(element, 'alfa_iD:', alfa_iD, 'alfa_iW', alfa_iW)

    return



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
        elementDicc['n-Octano']['value'].set(0.1)
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