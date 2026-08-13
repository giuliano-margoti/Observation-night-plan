from astroquery.jplhorizons import Horizons
import numpy as np
import time
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body
from sora import Body
from pathlib import Path
import pandas as pd

inicio = ['2027-02-16']
fim    = ['2027-02-23']

input_file = "objetos/OP2026B-003.txt"

df = pd.read_csv(input_file, header=None)

altura_max = 35.
mag_max = 22.
rot_max = 40.

site = 'OPD'
sulfix = 'Troianos_1m60'

time_init = " 21:00"
delta_time = 12

lista = np.array(df[0])#[20000588, 20000617, 20000624, 20000659, 20000884, 20001143]
print(lista)


    
rot = np.array([])
for spkid in lista:

    obj = Body(
            name=str(spkid),
            spkid=str(spkid),
            ephem=[f'bsp/{spkid}.bsp', 'bsp/de440.bsp'],
        )
    rot = np.append(rot, (obj.rotation.value))


for k in range(len(inicio)):
    
    masgs_list = np.array([])
    nome_obj = np.array([])
    
    for spkid in lista:
        obj = Horizons(
                    id=spkid,
                    location='500',
                    epochs={'start': inicio[k], 'stop': (Time(inicio[k])+1*u.min).value, 'step': f'{1}m'},
                    id_type='designation'
                )
        efem = obj.ephemerides(extra_precision=True)
        nome_obj = np.append(nome_obj, efem['targetname'][0])
        print(efem['targetname'][0])
        try:
            masgs_list = np.append(masgs_list,efem['V'][0])
        except:
            masgs_list = np.append(masgs_list,np.inf)
    
#       time.sleep(1)   

    if site == 'OPD':
    	location = EarthLocation(lat=-22.53502183332583*u.deg, lon=-45.58332688242411*u.deg, height=1849*u.m)

    elif site == 'SOAR':
    	location = EarthLocation(lat=-30.237806271993584*u.deg, lon=-70.73349396618212*u.deg, height=2674*u.m)
    contador_dia = 0
    times = Time(inicio[k]+' 21:00')
    while times < Time(fim[k]):
    
        obj = Horizons(
                        id='301',
                        location='500',
                        epochs={'start': times.value, 'stop': (times+1*u.min).value, 'step': f'{1}m'},
                    )
        
        efem = obj.ephemerides(extra_precision=True)
        
#        print(times)
    
        alt_list = {spkid: [] for spkid in lista}
        time_list = []
        time_local_list = []
        
        f_res = open("teste.txt", "w") 
        
        contador_hora = 0

        for contador_hora in range(delta_time):
        
            times = Time(inicio[k]+time_init) + contador_hora*u.h +contador_dia*u.day
            altaz_frame = AltAz(obstime=times, location=location)
            
            lst = times.sidereal_time('apparent', longitude=location.lon)
        
            lua = get_body("moon", times, location=location)
            lua_altaz = lua.transform_to(altaz_frame)
        
            local_time = times - 3*u.h
            f_res.write(f"\n\nLocal Time: {local_time.value[:16]} -- UTC: {times.value[:16]}\n")
        
            resultados = [] 
        
            time_local_list.append(local_time.value[11:16])
            time_list.append(times.value[11:16])

            for i, spkid in enumerate(lista):
        
                obj = Body(
                    name=spkid,
                    spkid=spkid,
                    ephem=[f'bsp/{spkid}.bsp', 'bsp/de440.bsp'],
                    database=None
                )
        
                pos = obj.get_position(times)
        
                coord = SkyCoord(ra=pos.ra, dec=pos.dec, frame='icrs')
                obj_altaz = coord.transform_to(altaz_frame)
        
                alt = obj_altaz.alt.deg
        
                if (alt > altura_max) & (masgs_list[i] < mag_max) & (rot[i] < rot_max):
        
                    alt_list[spkid].append(f"{alt:.1f}")
                    ra_hms = pos.ra.to_string(unit=u.hour, sep=' ', precision=2, pad=True)
                    dec_dms = pos.dec.to_string(unit=u.deg, sep=' ', precision=1, alwayssign=True, pad=True)
        
                    sep = obj_altaz.separation(lua_altaz)
                    
                    ha = (lst - coord.ra).wrap_at(24 * u.hourangle)
                    ha_formatado = ha.to_string(unit=u.hour, sep=':', precision=0, pad=True) 
        
                    resultados.append({
                        "nome": nome_obj[i],
                        "ra_dec": f"{ra_hms} {dec_dms}",
                        "alt": alt,
                        "az": obj_altaz.az.deg,
                        "ha": ha_formatado, 
                        "mag": masgs_list[i],
                        "moon_sep": sep.deg
                    })
                else:
                    
                    alt_list[spkid].append('')
        
            resultados_ordenados = sorted(resultados, key=lambda x: x["alt"], reverse=True)
        
            for obj in resultados_ordenados:
                f_res.write(
                    f"{obj['nome']:<35}\t| RA, DEC {obj['ra_dec']:<20}"
                    # --- NOVO: Inserido o AH na formatação do print ---
                    f"\t| EL, AZ, AH: {obj['alt']:<5.1f} {obj['az']:<5.1f} {obj['ha']:<8}"
                    f"\t| Mag: {obj['mag']:<5.1f}"
                    f"\t| moon_dis: {obj['moon_sep']:<5.1f}\n"
                )
            f_res.write("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
        
        
            print(times)
        
        f_res.close()
        
        
        obj_width = max(len('Objeto'), max(len(nome) for nome in nome_obj))
        time_widths = [max(len(t), 4) for t in time_list]
        
        
        dados = []
        for i, spkid in enumerate(lista):
        
            if all(val == '' for val in alt_list[spkid]):
                continue
        
            dados.append({
                "nome": nome_obj[i],
                "mag": masgs_list[i],
                "alturas": alt_list[spkid],
                "rot":rot[i]
            })

        
        dados_ordenados = sorted(dados, key=lambda x: x["mag"])
        
        obj_width = max(len('Objeto'), max(len(d["nome"]) for d in dados_ordenados))
        mag_width = 5
        time_widths = [max(len(t), 4) for t in time_list]
        
        with open("teste_tabela.txt", "w") as f:
        
            f.write("Tabela de observabilidade - Troianos - IAG\n\n\n")
        
            f.write(f"Iluninação da lua: {efem['lunar_illum'][0]:.1f} %\n\n\n")
        
            f.write("--------------------------------------------------------------------------------------------------------------------------------------------------------\n\n\n")
        
        
            header = f"{'Objeto':<{obj_width}} | {'mag':^{mag_width}} | {'rot':^{mag_width}} |"
            for i, t in enumerate(time_list):
                header += f" {t:^{time_widths[i]}} |"
            f.write(header + '\n')

            header = f"{'Objeto':<{obj_width}} | {'mag':^{mag_width}} | {'rot':^{mag_width}} |"
            for i, t in enumerate(time_local_list):
                header += f" {t:^{time_widths[i]}} |"
            f.write(header + '\n')
        
            separator = '-' * (obj_width + 1) + '+'
            separator += '-' * (mag_width + 2) + '+'
            separator += '-' * (mag_width + 2) + '+'
            for w in time_widths:
                separator += '-' * (w + 2) + '+'
            f.write(separator + '\n')
        
            for obj in dados_ordenados:
        
                line = f"{obj['nome']:<{obj_width}} | {obj['mag']:{mag_width}.1f} | {obj['rot']:{mag_width}.1f} |"
        
                for j, val in enumerate(obj["alturas"]):
                    if val:
                        line += f" {val:>{time_widths[j]}} |"
                    else:
                        line += f" {'':^{time_widths[j]}} |"
        
                f.write(line + '\n')
        
            f.write(separator + '\n\n\n')
            f.write("--------------------------------------------------------------------------------------------------------------------------------------------------------\n\n\n")
        
            f.write("Unidades: \nNome | RA (hh mm ss.ss), DEC (dd mm ss.s) | EL (dd.d) AZ (dd.d) AH (hh:mm:ss) | magnitude aparente | moon_dis (dd.d)")
        
        
        
        arquivo2 = Path("teste.txt")
        arquivo1 = Path("teste_tabela.txt")
        resultado = Path(f"plano_2026B/plano_da_noite_{sulfix}__{(times-1*u.day).value[:10]}.txt")
        
        conteudo = arquivo1.read_text(encoding='utf-8') + "\n" + arquivo2.read_text(encoding='utf-8')
        resultado.write_text(conteudo, encoding='utf-8')
        
        print(f"plano_2026B/plano_da_noite_{sulfix}__{(times-1*u.day).value[:10]}.txt")
        contador_dia+=1
