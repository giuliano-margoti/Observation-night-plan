from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import os


input_file = "filter.txt"
output_index = 'plano_da_noite_'

df = pd.read_csv(input_file, header=None)
tno_list = np.array(df[0])


t_inicio    = Time('2025-09-07 21:00')     
t_fim       = Time('2025-09-08 10:00')     
passo_min   = 60   
t0 = t_inicio.value[:10]
masgs_list = np.array([])

nome_obj = np.array([])

for nome_objeto in tno_list:
    print(nome_objeto)
    obj = Horizons(
                id=nome_objeto,
                location='500',
                epochs={'start': t_inicio.value, 'stop': (t_inicio+1*u.min).value, 'step': f'{passo_min}m'},
                id_type='designation'
            )
    efem = obj.ephemerides(extra_precision=True)
    nome_obj = np.append(nome_obj, efem['targetname'][0])
    #print(efem['targetname'][0])
    try:
        masgs_list = np.append(masgs_list,efem['V'][0])
    except:
        masgs_list = np.append(masgs_list,np.inf)
    
lista_ordenada = sorted(zip(tno_list, masgs_list, nome_obj), key=lambda item: item[1])

tno_list   = [item[0] for item in lista_ordenada]
masgs_list = [item[1] for item in lista_ordenada]
nome_obj   = [item[2] for item in lista_ordenada]



len_0 = 1#len(tno_list[0])

for nome_objeto in nome_obj:
    act = len(nome_objeto)
    if act > len_0:
        len_0 = act

cont = 0
for nome_objeto in nome_obj:
    act = len(nome_objeto)
    if act <= len_0:
        nome_obj[cont] =  nome_obj[cont]+' '*(len_0 - act + 1)
    print(f'|{nome_obj[cont]}| mag: {masgs_list[cont]}')
    cont+=1


separador_ = '------------------------------------------------------------------------------------------------------------------------'+'-'*len_0

code = '874'


# Filters
min_el = 35.


#with open(('plano_da_noite.txt', 'a+') as f_res:


##########################################################################################################################################################################################

alt_list = {}
nome_obj   = np.array([])
for nome_objeto in tno_list:
    alt_list[nome_objeto] = np.array([])


time_list = np.array([])

while t_inicio < t_fim:

    f_res = open("aux1.txt", "a+")

    local_time = t_inicio - 3*u.h
    time_list = np.append(time_list,local_time.value[11:16])
    print(f'\n\n{local_time}\n\n')
    
    f_res.write(f"\n\nLocal Time: {local_time.value[:16]} -- UTC: {t_inicio.value[:16]}\n")
    
    f_res.write(f"{separador_}\n")
    
    for nome_objeto in tno_list:
        
        efem = {}
        try:
            obj = Horizons(
                id=nome_objeto,
                location=code,
                epochs={'start': t_inicio.value, 'stop': (t_inicio+1*u.min).value, 'step': f'{passo_min}m'},
                id_type='designation'
            )
            efem = obj.ephemerides(extra_precision=True)
            act = len(efem['targetname'][0])

            if act <= len_0:
                nome_val =  efem['targetname'][0]+' '*(len_0 - act + 1)
            nome_obj = np.append(nome_obj, nome_val)
            mask =  (efem['EL'] > min_el)
            
            #time.sleep(0.5)
            
        except ValueError as e:

            if "Ambiguous target name; provide unique id" in str(e):
                print(f"Erro: {e}")
    
            else:

                print(f"Ocorreu um erro de valor: {e}")
    
        
        except Exception as e:

            print(f"Ocorreu um erro inesperado: {e}")
    
        try:
            mag_val = f"{efem['V'][mask][0]:.1f}"
        except:
            mag_val = np.inf
            
        try:
            
            ra_data = efem['RA'][mask].data
            dec_data = efem['DEC'][mask].data
            
            coords = SkyCoord(
                ra=ra_data * u.deg,
                dec=dec_data * u.deg,
                frame='icrs'
            )
            
            ra_hms = coords.ra.to_string(unit=u.hour, sep=' ', precision=2, pad=True)
            dec_dms = coords.dec.to_string(unit=u.deg, sep=' ', precision=1, alwayssign=True, pad=True)
            
            el_val = f"{efem['EL'][mask][0]:.1f}"
            moon_dis_val = f"{efem['lunar_elong'][mask][0]:.1f}"
            moon_lum_val = f"{efem['lunar_illum'][mask][0]:.1f}"
            
            
            
            f_res.write(f"{nome_val}\t| RA, DEC {ra_hms[0]} {dec_dms[0]}\t| EL: {el_val}\t| Mag: {mag_val}\t| moon_dis: {moon_dis_val}\t| moon_lum: {moon_lum_val}\n")
            alt_list[nome_objeto] = np.append(alt_list[nome_objeto], el_val)
        except:
            print(f"Sem dados para {efem['targetname'][0]}")
        
            alt_list[nome_objeto] = np.append(alt_list[nome_objeto], '')
    t_inicio = t_inicio+passo_min*u.min
    f_res.close()
    print('___________________________________________________________')
    




obj_width = max(len(obj) for obj in nome_obj)
time_widths = [max(len(t), len(max(alt_list[obj], key=len))) for t in time_list for obj in tno_list]


time_widths = [max(len(t), 4) for t in time_list] 
obj_width = max(len('Objeto'), max(len(obj) for obj in nome_obj))

with open("aux0.txt", "a+") as f_res:
    f_res.write(f"Tabela de observabilidade para a noite de {t0}\n\n\n{separador_}\n\nUnidades: RA, DEC hh mm ss.sss dd mm ss.s | EL: dd.d | Mag: mag aparente | moon_dis: dd.d | moon_lum: %\n\n\n")


    header = f"{'Objeto':<{obj_width}} |"
    for i, time in enumerate(time_list):
        header += f" {time:^{time_widths[i]}} |"
    f_res.write(header + '\n')


    separator = '-' * (obj_width + 1) + '+'
    for i in range(len(time_list)):
        separator += '-' * (time_widths[i] + 2) + '+'
    f_res.write(separator + '\n')
    contador = 0
    for nome_objeto in tno_list:
        
        line = f"{nome_obj[contador]:<{obj_width}} |"
        for i, time_val in enumerate(alt_list[nome_objeto]):

            if isinstance(time_val, str) and time_val.replace('.', '').isdigit():
                line += f" {time_val:>{time_widths[i]}} |"
            else:
                line += f" {time_val:^{time_widths[i]}} |"
        f_res.write(line + '\n')
        contador += 1
        

    f_res.write(separator + '\n')
        
    f_res.write(f'\n\n{separador_}\n\n\n')

        
arquivo1 = 'aux0.txt'
arquivo2 = 'aux1.txt'
arquivo_destino = f'{output_index}_{t0}.txt'

try:
    with open(arquivo_destino, 'w') as f_destino:
        
        with open(arquivo1, 'r') as f1:
            f_destino.write(f1.read())
            f_destino.write('\n')
        
        with open(arquivo2, 'r') as f2:
            f_destino.write(f2.read())
            
except Exception as e:
    print(f"Ocorreu um erro inesperado: {e}")
        
        

os.system('rm aux0.txt aux1.txt')
        
        
        
        
        
        
        
        
        
        
        
