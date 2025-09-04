from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

input_file = "rot.txt"

# Load the TNO IDs 
df = pd.read_csv(input_file, header=None)
tno_list = np.array(df[0])

code = '874'
t_inicio    = '2023-07-28 18:00'     
t_fim       = '2023-07-29 13:00'     
passo_min   = 10   

# Filters
min_el = 35.


for nome_objeto in tno_list:
    print(f'\n\n{nome_objeto}\n\n')
    efem = {}
    try:
        obj = Horizons(
            id=nome_objeto,
            location=code,
            epochs={'start': t_inicio, 'stop': t_fim, 'step': f'{passo_min}m'}
        )
        efem = obj.ephemerides(extra_precision=True)
        
        mask =  (efem['solar_presence'] != '*') & (efem['EL'] > min_el)
        
    except ValueError as e:

        if "Ambiguous target name; provide unique id" in str(e):
            print(f"Erro: {e}")

        else:

            print(f"Ocorreu um erro de valor: {e}")

    
    except Exception as e:

        print(f"Ocorreu um erro inesperado: {e}")

        
    try:
        datas_str = efem['datetime_str'][mask]
        datas_dt = [datetime.strptime(x, '%Y-%b-%d %H:%M') for x in datas_str]
        tempos = Time(datas_dt, scale='utc')
        
        datas_str = efem['datetime_str'][mask]

        datas_corrigidas = [
            (datetime.strptime(x, '%Y-%b-%d %H:%M') - timedelta(hours=3)).strftime('%Y-%b-%d %H:%M')
            for x in datas_str
        ]
    
        data = pd.DataFrame({
        'RA': pd.Series(efem['RA'][mask]).apply(lambda x: f"{x:.5f}"),
        'DEC': pd.Series(efem['DEC'][mask]).apply(lambda x: f"{x:.5f}"),
        'dst': pd.Series(efem['r'][mask]).apply(lambda x: f"{x:.3f}"),
        'time UTC': efem['datetime_str'][mask],
        'time LOCAL': datas_corrigidas,
        'mag': pd.Series(efem['V'][mask]).apply(lambda x: f"{x:.1f}"),
        'EL': pd.Series(efem['EL'][mask]).apply(lambda x: f"{x:.1f}"),
        'moon_dis': pd.Series(efem['lunar_elong'][mask]).apply(lambda x: f"{x:.1f}"),
        'moon_lum': pd.Series(efem['lunar_illum'][mask]).apply(lambda x: f"{x:.1f}")
    })
    

        ti = Time(datetime.strptime(datas_corrigidas[0], '%Y-%b-%d %H:%M'))
        tf = Time(datetime.strptime(datas_corrigidas[-1], '%Y-%b-%d %H:%M'))
        
        diff = (tf - ti).to(u.h)
        
        print(f'{nome_objeto} visivel por {diff}')
        
        data.to_csv('rot/'+nome_objeto.replace(' ', '')+'.txt',sep=' ', index=False, header=False)
    except:
        print(f'Sem dados para {nome_objeto}')
    
    print('___________________________________________________________')


