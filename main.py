# MAIN FILE
from tkinter import *
from functools import partial
from burbuja import temperaturaBurbuja
from flash_isotermico import start_flash_isotermico
from flash_adiabatico import start_flash_adiabatico
from termcolor import colored
from pprint import pprint
from copy import deepcopy
import math
import sys
import db

env_dev_flag = True

class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

def purgeDicc(zfDicc):
    new_dict = dict(zfDicc)
    for element in zfDicc:
        if(not zfDicc[element]['is_select'].get() == 1):
            del new_dict[element]
        else:
            new_dict[element]['Zi'] = float(zfDicc[element]['value'].get())
            del new_dict[element]['value']
            del new_dict[element]['is_select']
            if zfDicc[element]['is_xlf'].get() == 1:
                new_dict[element]['is_xlf'] = 1
                new_dict[element]['is_xhf'] = 0
            elif zfDicc[element]['is_xhf'].get() == 1:
                new_dict[element]['is_xhf'] = 1
                new_dict[element]['is_xlf'] = 0
            else:
                new_dict[element]['is_xlf'] = 0
                new_dict[element]['is_xhf'] = 0
    return new_dict

def calculate_TxxK(row_pos, columnaP):
    db_values = db.getElementValues(row_pos)
    A = db_values['A']
    B = db_values['B']
    C = db_values['C']
    temperatura = (B / (A - math.log(columnaP))) - C
    return temperatura


def start(f, zfDicc, p, tf, columnaP, xlk, xhk, R_OP):

    zfDicc = purgeDicc(zfDicc)
    f = float(f.get())
    p = float(p.get())
    tf = float(tf.get())
    columnaP = float(columnaP.get())
    xlk = float(xlk.get())
    xhk = float(xhk.get())
    R_OP = float(R_OP.get())

    # Determinar valores de pesado y ligero
    xlf = 0.0
    xlf_name = None
    xlf_list_pos = None
    ZLKF = None
    xhf = 0.0
    xhf_name = None
    xhf_list_pos = None
    ZHKF = None
    for element in zfDicc:
        zfDicc[element]['FxF'] = zfDicc[element]['Zi'] * f
        if zfDicc[element]['is_xlf'] == 1:
            xlf = zfDicc[element]['Zi']
            ZLKF = xlf
            xlf_name = element
            xlf_list_pos = zfDicc[element]['db_row']
        if zfDicc[element]['is_xhf'] == 1:
            xhf = zfDicc[element]['Zi']
            ZHKF = xhf
            xhf_name = element
            xhf_list_pos = zfDicc[element]['db_row']
    
    # Calculo de flujos de alimentacion de ligero y pesado
    flf = xlf * f
    fhf = xhf * f
	# Flujos de componentes en el domo y fondo
    fld = flf * xlk
    fhd = fhf * xhk
    fhw = fhf - fhd
    flw = flf - fld
    # Calculo de Ftd y Ftw
    Ftd = 0.0
    Ftw = 0.0
    i = 1
    for element in zfDicc:
        if i < xhf_list_pos:    
            fxf = zfDicc[element]['FxF']
            Ftd += fxf 
        elif i >= xhf_list_pos:
            fxf = zfDicc[element]['FxF']
            Ftw += fxf
        i += 1
    
    Ftd += fld + fhd
    Ftw += flw + fhw
    

    # Valores de arranque para Domo
    i = 1
    for element in zfDicc:
        if i < xhf_list_pos and i != xlf_list_pos:
            xi = zfDicc[element]['FxF'] / Ftd
            zfDicc[element]['Xi'] = xi
        elif i == xlf_list_pos:
            xi = fld / Ftd
            zfDicc[element]['Xi'] = xi
        elif i == xhf_list_pos:
            xi = fhd / Ftd
            zfDicc[element]['Xi'] = xi
        else:
            zfDicc[element]['Xi'] = 0
        i += 1
    
    print('Arranque para D')
    for element in zfDicc:
        print(element, ':', zfDicc[element]['Xi'])
    # Temperatura Domo
    TDLK = calculate_TxxK(xlf_list_pos, columnaP)
    # Temperatura W
    TWHK = calculate_TxxK(xhf_list_pos, columnaP)
    print('Ligero:', xlf_name, 'Pesado:', xhf_name)
    print('TDLK:', TDLK, 'TWHK:', TWHK)
    print('Ftd:', Ftd, 'Ftw:', Ftw)
    # Copia de diccionario
    zfDiccV2 = deepcopy(zfDicc)
    # Valores de arranque para fondo
    i = 1
    for element in zfDiccV2:
        if i > xhf_list_pos:
            xiw = zfDiccV2[element]['FxF'] / Ftw
            zfDiccV2[element]['Xi'] = xiw
        elif i == xhf_list_pos:
            xiw = fhw / Ftw
            zfDiccV2[element]['Xi'] = xiw
        elif i == xlf_list_pos:
            xiw = flw / Ftw
            zfDiccV2[element]['Xi'] = xiw
        else:
            zfDiccV2[element]['Xi'] = 0
        i += 1
    print('Arranque para W')
    for element in zfDiccV2:
        print(element, ':', zfDiccV2[element]['Xi'])

    megaiterator = 0
    dictCompare = None
    while True:
        print(" --------- TEMP BURBUJA Domo ------------")

        it_count = 0
        result_d_burbuja = temperaturaBurbuja(columnaP, zfDicc, TDLK)
        while not result_d_burbuja['status']:
            it_count+=1
            print('Iteracion en D:', it_count)
            result_d_burbuja = temperaturaBurbuja(columnaP, zfDicc, result_d_burbuja['Td'])
        
        Td_obtenida = result_d_burbuja['Td']
        kibs_d_obtenida = result_d_burbuja['kib_calculadas']
       
        
        print(" --------- TEMP BURBUJA W ------------")
        it_count = 0
        result_w_burbuja = temperaturaBurbuja(columnaP, zfDiccV2, TWHK)
        while not result_w_burbuja['status']:
            it_count+=1
            print('Iteracion en W:', it_count)
            result_w_burbuja = temperaturaBurbuja(columnaP, zfDiccV2, result_w_burbuja['Td'])
        
        Tw_obtenida = result_w_burbuja['Td']
        kibs_w_obtenida = result_w_burbuja['kib_calculadas']
        
        print('TW Obtenida:', Tw_obtenida)
        print('TD Obtenida:', Td_obtenida)

        alfa_iD = 0.0
        alfa_iW = 0.0
        print('Kibs D:')
        pprint(kibs_d_obtenida)
        print('Kibs W:')
        pprint(kibs_w_obtenida)
        KidH = kibs_d_obtenida[xhf_name]
        KiwH = kibs_w_obtenida[xhf_name]
        for element in zfDicc:
            alfa_iD = kibs_d_obtenida[element] / KidH
            alfa_iW = kibs_w_obtenida[element] / KiwH
            zfDicc[element]['alfa_iD'] = alfa_iD
            zfDicc[element]['alfa_iW'] = alfa_iW
            print(element, 'alfa_iD:', alfa_iD, 'alfa_iW', alfa_iW)

        alfaLKD = zfDicc[xlf_name]['alfa_iD']
        alfaLKW = zfDicc[xlf_name]['alfa_iW']
        alfaL_prom = math.sqrt(alfaLKW * alfaLKW)

        Nmin = math.log((fld / fhd) * (fhw / flw)) / math.log(alfaL_prom)

        for element in zfDicc:
            alfaID_prom = math.sqrt(kibs_d_obtenida[element] * zfDicc[element]['alfa_iD'])
            alfaIW_prom = math.sqrt(kibs_w_obtenida[element] * zfDicc[element]['alfa_iW'])
            zfDicc[element]['alfaID_prom'] = alfaID_prom
            zfDicc[element]['alfaIW_prom'] = alfaIW_prom

        D = 0.0
        W = 0.0
        for element in zfDicc:
            if element != xlf_name and element != xhf_name:
                F = zfDicc[element]['FxF']
                alfaID_prom = zfDicc[element]['alfaID_prom']
                alfaIW_prom = zfDicc[element]['alfaIW_prom'] 
                biNK = F / (1 + ((fhd / fhw) * alfaIW_prom**Nmin))
                diNK = (F * (fhd / fhw) * alfaID_prom**Nmin) / (1 + ((fhd / fhw) * alfaID_prom**Nmin))
                zfDicc[element]['biNK'] = biNK
                zfDicc[element]['diNK'] = diNK
                if diNK > 0:
                    D += diNK
                if biNK > 0:
                    W += biNK

        D += fld + fhd
        W += fhw + flw

        if megaiterator > 0:
            if comparate(zfDicc, dictCompare, xlf_name, xhf_name):
                break
        else:
            dictCompare = dict(zfDicc)
            TDLK = Td_obtenida
            TWHK = Tw_obtenida
        
        for element in zfDicc:
            if element != xlf_name and element != xhf_name:
                zfDicc[element]['Xi'] = zfDicc[element]['diNK'] / D 
                zfDiccV2[element]['Xi'] = zfDicc[element]['biNK'] / W
            elif element == xlf_name:
                zfDicc[element]['Xi'] = fld / D # D
                zfDiccV2[element]['Xi'] = flw / W # W
            elif element == xhf_name:
                zfDicc[element]['Xi'] = fhd / D # D
                zfDiccV2[element]['Xi'] = fhw / W # W

        megaiterator+=1
    print(colored('______________________COMPLETE ITERATION LOOP AT: ' + str(megaiterator), 'red'))

    print(" --------- FLASH ISOTERMICO------------")
    result_flash_iso = start_flash_isotermico(tf, p, zfDicc)
    print(" --------- FLASH ADIABATICO------------")
    result_flash_adia = start_flash_adiabatico(tf, columnaP, zfDicc)
    Hv_A = result_flash_adia['Hv_A']
    HL_A = result_flash_adia['HL_A']
    HF = result_flash_iso['HF']
    psi_A = result_flash_adia['psi_A']
    psi_F = result_flash_iso['psi']
    fTvAd = (psi_A * Hv_A) + ((1 - psi_F)*HL_A) - HF
    fTvAd = fTvAd / 1000
    while fTvAd >= -0.0001 and fTvAd <= 0.0001:
        break
        TvAd = result_flash_adia['TvAd']
        if TvAd > 0:
            TvAd -= 0.1
        elif TvAd < 0:
            TvAd += 0.1
        result_flash_adia = start_flash_adiabatico(TvAd, columnaP, zfDicc)
        fTvAd = (psi_A * Hv_A) + ((1 - psi_F)*HL_A) - HF
        fTvAd = fTvAd / 1000

    print(" --------- UNDERWOOD ------------")
    q = (Hv_A - HF) / (Hv_A - HL_A)
    for element in zfDicc:
        zfDicc[element]['alfa_i'] = zfDicc[element]['Kiad'] / zfDicc[xhf_name]['Kiad']

    Rmin = 0.0
    if columnaP < 4:
        Rmin = 0.0003*p**4 - 0.0078*p**3 + 0.0865*p**2 - 0.4219*p + 1.5406
    elif columnaP >= 4 and columnaP < 8:
        Rmin = 0.0003*p**4 - 0.0078*p**3 + 0.0865*p**2 - 0.4219*p + 1.7406
    elif columnaP >= 8:
        Rmin = 0.0003*p**4 - 0.0078*p**3 + 0.0865*p**2 - 0.4219*p + 1.9406

    print(" --------- Gillilaand ------------")
    BIGX = (R_OP - Rmin) / (R_OP + 1)
    N = ((1 - math.exp((1 + (54.4)*BIGX) / (11 + (117.2*BIGX)))) * (BIGX - 1 / BIGX**0.5)) + Nmin
    N = N / (1 - (1 - math.exp((1 + (54.4 * BIGX)) /  (11 + (117.2*BIGX))))) * (BIGX - 1 / BIGX**0.5)

    print(colored(' --------- RESULTADO FINAL ------------', 'cyan'))
    NR_NS_rel = ((ZHKF/ZLKF) * (W/D) * ((xlk * W) / (xhk * D)))**0.206
    print(NR_NS_rel)
    # Zona de Retificación:
    NR = (NR_NS_rel / (1 + NR_NS_rel)) * N
    NS = NR - N
    print('NR_NS_rel:', NR_NS_rel, 'NR:', NR, 'NS:', NS)
    return NR_NS_rel

def comparate(dict1, dict2, xlf_name, xhf_name):
    for element in dict1:
        if element != xlf_name and element != xhf_name:
            biNK1 = dict1[element]['biNK']
            diNK1 = dict1[element]['diNK']
            biNK2 = dict2[element]['biNK']
            diNK2 = dict2[element]['diNK']
            if abs(biNK1 - biNK2) > 0.1:
                return False
            if abs(diNK1 - diNK2) > 0.1:
                return False
    return True

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
    R_OP = StringVar()
    xhkEntry = Entry(root, textvariable=R_OP).grid(row=rowCounter, column=1)

    if(env_dev_flag):
        f.set(100)
        for element in elementDicc:
            elementDicc[element]['is_select'].set(1)
            # if element == 'n-Hexano':
            #     print('----->   Set Etano to XLF',)
            #     elementDicc[element]['is_xlf'].set(1)
            # if element == 'Acetato de Etilo':
            #     print('----->   Set Cumeno to XHF',)
            #     elementDicc[element]['is_xhf'].set(1)

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
        tf.set(120)
        columnaP.set(4)
        xlk.set(0.97)
        xhk.set(0.03)
        R_OP.set(0.4)
    startsSimulation = partial(start, f, elementDicc, p, tf, columnaP, xlk, xhk, R_OP)

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