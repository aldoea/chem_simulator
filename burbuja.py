import numpy as np
import math
import db

def temperaturaBurbuja(p, zfDicc, tf):
    Td_obtenida = 0.0
    yis_calculadas = list()
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
        Kib = calculate_Kib(db_values, zfDicc, td, p)
        while True:
            # Calcular yi_calculada
            yi_calculada=Kib*xi
            #print('yi calculada = ', Kib,'*',xi, '=', yi_calculada)
            normalizado = abs(yi_supuesta - yi_calculada)
            #print('abs(yi_supuesta - yi_calculada) => ', yi_supuesta, '-', yi_calculada, '=', normalizado)
            if(normalizado  < 0.001):
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
                Kib = calculate_Kib(db_values, zfDicc, td, p)
                print('Nuevo Kib:', Kib)

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

def calculate_Kib(db_values, zfDicc, tf, p):
    Tc = db_values['Tc']
    Pc = db_values['Pc']
    Ai = calculate_Ai(Pc, tf, Tc)
    Ay = calculate_Ay(zfDicc, tf, p)
    Ax = calculate_Ax(zfDicc, tf)
    Bi = calculate_Bi(Pc, tf, Tc)
    By = calculate_By(zfDicc, tf, p)
    Bx = calculate_Bx(zfDicc, tf)

    A = math.sqrt(Ax * Ay)
    B = math.sqrt(Bx * By)

    Z = getRoots(A, B, p)
    Zl = Z['Zl']
    Zv = Z['Zv']
    FIv = coeficiente_de_fugacidad(Zv, p, A, B, Ai, Bi)
    FIl = coeficiente_de_fugacidad(Zl, p, A, B, Ai, Bi)
    #print('FIv:', FIv)
    #print('FIl:', FIl)
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
        xi = float(zfDicc[element]['value'].get())
        sumatoria += Bi * xi
    return sumatoria

def calculate_Ax(zfDicc, tf):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Ai = calculate_Ai(Pc, tf, Tc)
        xi = float(zfDicc[element]['value'].get())
        sumatoria += Ai * xi
    return sumatoria

def calculate_Ay(zfDicc, tf, p):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Ai = calculate_Ai(Pc, tf, Tc)
        xi = float(zfDicc[element]['value'].get())
        ki = calculate_Ki(p, tf, db_values)
        yi_supuesta = calculate_yi_supuesta(ki, xi)
        sumatoria += Ai * yi_supuesta
    return sumatoria

def calculate_By(zfDicc, tf, p):
    sumatoria = 0.0
    for element in zfDicc:
        db_values = db.getElementValues(zfDicc[element]['db_row'])
        Tc = db_values['Tc']
        Pc = db_values['Pc']
        Bi = calculate_Bi(Pc, tf, Tc)
        xi = float(zfDicc[element]['value'].get())
        ki = calculate_Ki(p, tf, db_values)
        yi_supuesta = calculate_yi_supuesta(ki, xi)
        sumatoria += Bi * yi_supuesta
    return sumatoria

def getRoots(A, B, p):
    coef_c = (B * p) * ( (A**2/B) - (B*p) - 1)
    coef_d = -1 * ((A**2/B) * ((B*p)**2))
    roots = np.roots([1, -1, coef_c, coef_d])
    # for i in range(0, len(roots)):
    #     roots[i] = float(roots[i])
    roots.sort()
    #print('Raices:', roots)
    Z = dict()
    Z['Zl'] = roots[2]
    Z['Zv'] = roots[1]
    #print('Zv:',Z['Zv'])
    #print('Zl:',Z['Zl'])
    return Z

def coeficiente_de_fugacidad(Z, p, A, B, Ai, Bi):
    first_eq = (Z - 1) * (Bi / B)
    second_eq = math.log(Z - (B * p))
    subeq3_1 = (A**2 / B) * ((2*Ai / A) - (Bi / B))
    subeq3_2 = math.log(1 + ((B * p) / Z))
    third_eq = subeq3_1 * subeq3_2
    final_eq = first_eq - second_eq - third_eq
    return math.exp(final_eq)