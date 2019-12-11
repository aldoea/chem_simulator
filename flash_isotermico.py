import db
import random
import math
import burbuja
import sys
from termcolor import colored, cprint

R = 8.314*10**-5
H0_iv = [0.05, 0.06, 0.08, 0.9, 0.4, 0.02, 0.03, 0.05, 0.06, 0.09]
TEMP_CORRECTOR = 273.15
def scitter():
    input()
    sys.exit()

def start_flash_isotermico(tf, p, zfDicc):
    # Set TD
    td = tf
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Zif = zfDicc[element]['Zi']
        Kils = burbuja.calculate_Ki(p, td, db_values)
        zfDicc[element]['Kils'] = Kils
        zfDicc[element]['Xils'] = Zif
        print(element,'Zi:', Zif, 'Kils:', Kils)
        zfDicc[element]['Yils'] = Zif * Kils
        zfDicc[element]['is_normalized'] = False

    psi = 0.5
    normalized_counter = 0
    psi_kplus1 = recalculate_psi(zfDicc, psi)
    while True:
        psi_abs = abs((psi_kplus1 - psi) / psi_kplus1)
        print('psi_abs:', psi_abs)
        if psi_abs < 0.001:
            print('©- psi_abs pass the error prove, Psi:', psi, 'psi^k+1:', psi_kplus1)
            #scitter()
            for element in zfDicc:
                print(colored(element, 'yellow'))
                if not zfDicc[element]['is_normalized']:
                    db_values = db.getElementValues(zfDicc[element]['db_row'])
                    Zif = zfDicc[element]['Zi']
                    Kils = zfDicc[element]['Kils']
                    xils_supuesta = zfDicc[element]['Xils']
                    yils_supuesta = zfDicc[element]['Yils']
                    xil_calculada = Zif / (1 + (psi_kplus1 * (Kils - 1)))
                    yil_calculada = (Zif * Kils) / (1 + (psi_kplus1 * (Kils - 1)))
                    norm_xi = abs(xils_supuesta - xil_calculada)
                    norm_yi = abs(yils_supuesta - yil_calculada)
                    if norm_xi < 0.001 and norm_yi < 0.001:
                        zfDicc[element]['Xil'] = xil_calculada
                        zfDicc[element]['Yil'] = yil_calculada
                        zfDicc[element]['is_normalized'] = True
                        zfDicc[element]['psi'] = psi_kplus1
                        print(colored('<-- ' + element + 'is normalized !!!', 'cyan'))
                        normalized_counter+=1
                        break
                    else:
                        print(colored('-- NOT NORMALIZED ' + element, 'red'))
                        print('-> Recalculating Kils for', element, ':', Kils)
                        zfDicc[element]['Xils'] = xil_calculada
                        zfDicc[element]['Yils'] = yil_calculada
                        valores_de_arranque = calculate_valores_de_arranque(db_values, zfDicc, p, td)
                        Kils = valores_de_arranque['Kils']
                        zfDicc[element]['Kils'] = Kils
                        print('-> Now Kils for', element, ':', Kils)
                        print('-> psi then:', psi, 'psi^k+1 then:', psi_kplus1)
                        psi = 0.5
                        psi_kplus1 = recalculate_psi(zfDicc, psi)
                        print('-> psi now:', psi, 'psi^k+1 now:', psi_kplus1)
                        # normalized_counter = 0
                        break
        else:
            print('-> psi then:', psi, 'psi^k+1 then :', psi_kplus1)
            psi = psi_kplus1
            psi_kplus1 = recalculate_psi(zfDicc, psi)
            print('-> psi now:', psi, 'psi^k+1 now:', psi_kplus1)
        if normalized_counter == len(zfDicc):
            break
    HF = calculate_HF(zfDicc, td, p)
    psi_avg = 0.0
    for element in zfDicc:
        psi_avg += zfDicc[element]['psi']
    psi_avg = abs(psi_avg / len(zfDicc))
    print('PSI adecuado:', psi_avg)
    resultado = dict()
    resultado['HF'] = HF
    resultado['psi'] = psi_avg
    return resultado

def recalculate_psi(zfDicc, psi):
    sum1 = 0.0
    sum2 = 0.0
    for element in zfDicc:
        Zif = zfDicc[element]['Zi']
        Kils = zfDicc[element]['Kils']
        tmp1 = (Zif * (1 - Kils)) / (1 + (psi * (1 - Kils)))
        tmp2 = (Zif * (1 - Kils)**2) / (1 + (psi * (1 - Kils)))**2
        sum1 += tmp1
        sum2 += tmp2
    new_psi = psi - (sum1 / sum2)
    return new_psi

def calculate_HF(zfDicc, tf, p):
    Ax = 0.0
    Ay = 0.0
    Bx = 0.0
    By = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Ai = burbuja.calculate_Ai(Pc, tf, Tc)
        Bi = burbuja.calculate_Bi(Pc, tf, Tc)
        Xil = zfDicc[element]['Xil']
        Yil = zfDicc[element]['Yil']
        Bx += Bi * Xil
        By += Bi * Yil
        Ax += Ai * Xil
        Ay += Ai * Yil

    A = math.sqrt(Ax * Ay)
    B = math.sqrt(Bx * By)

    Z = burbuja.getRoots(A, B, p)
    Zl = Z['Zl']
    Zv = Z['Zv']
    Hv = calculate_H(Zv, tf, p, 'Yil', A, B, zfDicc)
    HL = calculate_H(Zl, tf, p, 'Xil', A, B, zfDicc)
    print('Hv:', Hv)
    print('HL:', HL)
    print('HF:', Hv + HL)
    return Hv + HL

def calculate_H(Zx, tf, p, key, A, B, zfDicc):
    sumatoria = 0.0
    iterator = 0
    tf = (tf + TEMP_CORRECTOR) * 1.8
    for element in zfDicc:
        ii = zfDicc[element][key]
        tmp = ii * H0_iv[iterator]
        iterator+=1
        sumatoria += tmp
    second_eq = (R*tf) * (Zx - 1 - ((3*(A**2)) / (2*B)) * math.log(1 + ((B*p) / Zx)))
    return sumatoria + second_eq

def calculate_valores_de_arranque(db_values, zfDicc , p, tf):
    Tc = db_values['Tc']
    Tc = (Tc + TEMP_CORRECTOR) * 1.8
    Pc = db_values['Pc']
    Ax = 0.0
    Ay = 0.0
    Bx = 0.0
    By = 0.0
    tf = (tf + TEMP_CORRECTOR) * 1.8
    Aix = calculate_Ai(Pc, tf, Tc)
    Bix = calculate_Bi(Pc, tf, Tc)
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc_e = db_values['Tc']
        Tc_e = (Tc_e + TEMP_CORRECTOR) * 1.8
        Pc = db_values['Pc']
        Ai = calculate_Ai(Pc, tf, Tc_e)
        Bi = calculate_Bi(Pc, tf, Tc_e)
        Xil = zfDicc[element]['Xils']
        Yil = zfDicc[element]['Yils']
        Bx += (Bi * Xil)
        By += (Bi * Yil)
        Ax += (Ai * Xil)
        Ay += (Ai * Yil)

    A = math.sqrt(Ax * Ay)
    B = math.sqrt(Bx * By)

    valores_de_arranque = burbuja.getRoots(A, B, p)
    Zl = valores_de_arranque['Zl']
    Zv = valores_de_arranque['Zv']
    FIv = burbuja.coeficiente_de_fugacidad(Zv, p, A, B, Aix, Bix)
    FIl = burbuja.coeficiente_de_fugacidad(Zl, p, A, B, Aix, Bix)
    valores_de_arranque['FIv:'] = FIv
    valores_de_arranque['FIl:'] = FIl
    Kils = FIl/FIv
    valores_de_arranque['Kils'] = Kils
    return valores_de_arranque


def calculate_Ai(Pc, tf, Tc):
    return (6.20449/(Pc*((tf/Tc)**2.5)))**(0.5)

def calculate_Bi(Pc, tf, Tc):
    return 1.2574/(Pc*(tf/Tc))      
