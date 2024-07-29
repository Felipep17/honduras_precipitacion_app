# Standard library imports
import os
import time
import ssl
import shutil
import glob
import itertools
from datetime import datetime, timedelta
import io
import json
import math
from concurrent.futures import ThreadPoolExecutor
from tkinter import messagebox, filedialog, ttk  # tk and ttk for GUI development
import tkinter as tk
# Data processing and numerical analysis
import numpy as np
import pandas as pd
# Xarray and related libraries
import xarray as xr
import rioxarray as rxr
import xcast as xc
# Geospatial data handling
import geopandas as gpd
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# Plotting
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
# Statistical analysis
from scipy.stats import kendalltau, pearsonr
from scipy.linalg import svd
# Progress bar for loops
from tqdm import tqdm
# Requests for web access
import requests
def get_gcm_url(start_time, start_year, last_year, start_step, last_step, lat_1, lat_2, lon_1, lon_2, model):
    """
    Esta función genera una URL para descargar datos de precipitación de diferentes modelos GCM.

    Parámetros:
    - start_time (str): La hora de inicio en formato 'HHMM'.
    - start_year (int): El primer año del período de datos.
    - last_year (int): El último año del período de datos.
    - start_step (int): El paso de inicio (en días).
    - last_step (int): El paso final (en días).
    - lat_1 (float): La latitud inicial del área de interés.
    - lat_2 (float): La latitud final del área de interés.
    - lon_1 (float): La longitud inicial del área de interés.
    - lon_2 (float): La longitud final del área de interés.
    - model (str): El nombre del modelo GCM a utilizar. Debe ser una de las siguientes opciones:
        'NCEP-CFSv2', 'Can-1', 'Can-2', 'CCSM4', 'GFDL', 'NASA', 'ECMWF', 'Meteo-France', 
        'Glosea', 'CMCC', 'DWD', 'JMA'.
    
    Retorna:
    - str: La URL generada para el modelo especificado.
    
    Lanza:
    - ValueError: Si el modelo especificado no es válido.
    """
    
    # Diccionario que contiene las URLs base para cada modelo
    dic = {
        'NCEP-CFSv2': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NCEP-CFSv2/.HINDCAST/.PENTAD_SAMPLES_FULL/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-" + str(last_year) + ")/VALUES/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/M/1/28/RANGE/%5BM%5Daverage/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'Can-1': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC3/.HINDCAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/SOURCES/.Models/.NMME/.CanSIPS-IC3/.FORECAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/M/1/20/RANGE/%5BM%5Daverage/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'Can-2': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.CanSIPS-IC3/.CanCM4i-IC3/.HINDCAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/SOURCES/.Models/.NMME/.CanSIPS-IC3/.CanCM4i-IC3/.FORECAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/M/1/10/RANGE/%5BM%5Daverage/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'CCSM4': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.COLA-RSMAS-CCSM4/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-" + str(last_year) + ")/VALUES/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/M/1/10/RANGE/%5BM%5Daverage/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'GFDL': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.GFDL-SPEAR/.HINDCAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2019)/VALUES/M/1/15/RANGE/%5BM%5Daverage/SOURCES/.Models/.NMME/.GFDL-SPEAR/.FORECAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%202020-" + str(last_year) + ")/VALUES/M/1/30/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'NASA': "https://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.NASA-GEOSS2S/.HINDCAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2016)/VALUES/M/1/4/RANGE/%5BM%5Daverage/SOURCES/.Models/.NMME/.NASA-GEOSS2S/.FORECAST/.MONTHLY/.prec/S/(0000%201%20" + start_time + "%202017-" + str(last_year) + ")/VALUES/M/1/10/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'ECMWF': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.ECMWF/.SEAS51/.hindcast/.prcp/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2021)/VALUES/M/1/25/RANGE/[M]average/SOURCES/.EU/.Copernicus/.CDS/.C3S/.ECMWF/.SEAS51/.forecast/.prcp/S/(0000%201%20" + start_time + "%202022-" + str(last_year) + ")/VALUES/M/1/51/RANGE/[M]average/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/[L]//keepgrids/average/%28mm/day%29/unitconvert/c:/90//units//days/def/:c/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'Meteo-France': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.Meteo_France/.System8/.hindcast/.prcp/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2021)/VALUES/M/1/25/RANGE/%5BM%5Daverage/SOURCES/.EU/.Copernicus/.CDS/.C3S/.Meteo_France/.System8/.forecast/.prcp/S/(0000%201%20" + start_time + "%202022-" + str(last_year) + ")/VALUES/M/1/51/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/%28mm/day%29/unitconvert/c:/90//units//days/def/:c/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'Glosea': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.UKMO/.GloSea6-GC2/.System600/.hindcast/.prcp/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/M/1/28/RANGE/%5BM%5Daverage/SOURCES/.EU/.Copernicus/.CDS/.C3S/.UKMO/.GloSea6-GC2/.System600/.forecast/.prcp/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/M/1/50/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/(mm/day)/unitconvert/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'CMCC': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.CMCC/.SPSv3p5/.hindcast/.prcp/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/M/1/40/RANGE/%5BM%5Daverage/SOURCES/.EU/.Copernicus/.CDS/.C3S/.CMCC/.SPSv3p5/.forecast/.prcp/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/M/1/50/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/(mm/day)/unitconvert/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'DWD': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.DWD/.GCFS2p1/.hindcast/.prcp/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/M/1/30/RANGE/%5BM%5Daverage/SOURCES/.EU/.Copernicus/.CDS/.C3S/.DWD/.GCFS2p1/.forecast/.prcp/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/M/1/50/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/(mm/day)/unitconvert/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc",
        'JMA': "https://iridl.ldeo.columbia.edu/SOURCES/.EU/.Copernicus/.CDS/.C3S/.JMA/.CPS2/.hindcast/.prec/S/(0000%201%20" + start_time + "%20" + str(start_year) + "-2020)/VALUES/M/1/10/RANGE/%5BM%5Daverage/SOURCES/.EU/.Copernicus/.CDS/.C3S/.JMA/.CPS2/.forecast/.prec/S/(0000%201%20" + start_time + "%202021-" + str(last_year) + ")/VALUES/M/1/91/RANGE/%5BM%5Daverage/appendstream/L/" + str(start_step) + "/" + str(last_step) + "/RANGE/%5BL%5D//keepgrids/average/(mm/day)/unitconvert/c%3A/90//units//days/def/%3Ac/mul/Y/(" + str(lat_1) + ")/(" + str(lat_2) + ")/RANGEEDGES/X/(" + str(lon_1) + ")/(" + str(lon_2) + ")/RANGEEDGES/-999/setmissing_value/data.nc"
    }
    
    # Verificar si el modelo especificado está en el diccionario
    if model not in dic:
        raise ValueError("Modelo no válido. Debe ser uno de: 'NCEP-CFSv2', 'Can-1', 'Can-2', 'CCSM4', 'GFDL', 'NASA', 'ECMWF', 'Meteo-France', 'Glosea', 'CMCC', 'DWD', 'JMA'.")

    # Retornar la URL correspondiente al modelo especificado
    return dic[model]
def download_nc_file(url, path, timeout=300):
    """
    Descarga un archivo NetCDF desde una URL y lo guarda en la ruta especificada.
    """
    username = 'diegoagudelo30@gmail.com'
    password = 'diego2020'
    login_url = 'https://iridl.ldeo.columbia.edu/auth/login/local/submit/login'

    try:
        session = requests.Session()
        payload = {'email': username, 'password': password, 'redirect': 'https://iridl.ldeo.columbia.edu/auth'}
        login_response = session.post(login_url, data=payload)
        if login_response.status_code != 200 or 'login' in login_response.url:
            print(f"Error al iniciar sesión: {login_response.status_code}")
            return False

        print("Inicio de sesión exitoso")
        time.sleep(6)

        response = session.get(url, stream=True, timeout=timeout, verify=True)
        response.raise_for_status()

        expected_size = int(response.headers.get('Content-Length', 0))
        downloaded_size = 0
        progress_bar = tqdm(total=expected_size, unit='B', unit_scale=True, desc=path.split('/')[-1], miniters=1)

        with open(path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    downloaded_size += len(chunk)
                    f.write(chunk)
                    progress_bar.update(len(chunk))
        progress_bar.close()

        if downloaded_size == expected_size:
            print(f"Archivo descargado exitosamente: {path}")
            return True
        else:
            print(f"El tamaño del archivo descargado ({downloaded_size} bytes) no coincide con el tamaño esperado ({expected_size} bytes).")
            if os.path.exists(path):
                os.remove(path)
    except (requests.exceptions.RequestException, IOError) as e:
        print(f"Error en la descarga: {e}")
        if os.path.exists(path):
            os.remove(path)

    return False

def generate_chirps_url(start_day, start_month, start_year, end_day, end_month, end_year, lat_min, lat_max, lon_min, lon_max):
    """
    Genera una URL para descargar datos de precipitación diaria mejorada de CHIRPS.
    """
    base_url = "https://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.daily-improved/.global/.0p25/.prcp/T"
    start_date = f"({start_month}%{start_day}{start_year})"
    end_date = f"({end_month}%{end_day}{end_year})"
    lat_range = f"/Y/({lat_min})/({lat_max})/RANGEEDGES"
    lon_range = f"/X/({lon_min})/({lon_max})/RANGEEDGES"
    url = f"{base_url}/{start_date}/{end_date}/RANGE{lat_range}{lon_range}/data.nc"
    return url

def calculate_winter_precip(path, months, start_year, end_year):
    """
    Calcula la precipitación total acumulada para los meses de invierno especificados
    entre los años de inicio y fin proporcionados.
    """
    try:
        precip_chirps_daily = xr.open_dataset(path, decode_times=False).prcp.expand_dims({'M': [0]}, axis=1)
        precip_chirps_daily['T'] = pd.to_datetime(precip_chirps_daily['T'], unit='D', origin='julian')
        results = []

        for year in range(start_year, end_year + 1):
            start_month = months[0]
            end_month = months[-1]

            if start_month <= end_month:
                start_date = f'{year}-{int(start_month):02d}-01'
                end_date = f'{year}-{int(end_month):02d}-{pd.Timestamp(f"{year}-{int(end_month):02d}-01").days_in_month}'
                save_year = year
            else:
                start_date = f'{year}-{int(start_month):02d}-01'
                end_date = f'{year + 1}-{int(end_month):02d}-{pd.Timestamp(f"{year + 1}-{int(end_month):02d}-01").days_in_month}'
                save_year = year + 1

            filtered_precip_chirps_daily = precip_chirps_daily.sel(T=slice(start_date, end_date))
            available_dates = pd.to_datetime(filtered_precip_chirps_daily['T'].values)
            months_present = set(available_dates.month)
            required_months = set(pd.date_range(start=start_date, end=end_date, freq='MS').month)

            if required_months.issubset(months_present):
                total_precipitation = filtered_precip_chirps_daily.sum(dim='T')
                total_precipitation = total_precipitation.expand_dims(year=[save_year])
                results.append(total_precipitation)
            else:
                print(f"No se encontraron todos los meses requeridos para el año {year}.")

        if results:
            winter_precip = xr.concat(results, dim='year')
            return winter_precip
        else:
            print('No se encontraron datos completos para el rango de fechas especificado.')
            return None
    except FileNotFoundError:
        print(f"El archivo en la ruta {path} no fue encontrado.")
        return None
    except Exception as e:
        print(f"Ocurrió un error durante el cálculo de precipitación: {e}")
        return None

def convert_lat_lon(value):
    if 'N' in value:
        return float(value.replace('N', ''))
    elif 'S' in value:
        return -float(value.replace('S', ''))
    elif 'E' in value:
        return float(value.replace('E', ''))
    elif 'W' in value:
        return -float(value.replace('W', ''))
    return float(value)

def convert_month(month_str):
    months = {
        'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4,
        'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8,
        'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
    }
    return months.get(month_str, 1)

def kendall_tau(a, b):
    """
    Calcula el coeficiente de correlación de Kendall Tau entre dos arrays.
    """
    try:
        if len(a) != len(b):
            raise ValueError("Los arrays de entrada deben tener la misma longitud.")
        tau, _ = kendalltau(a, b)
        return tau
    except ValueError as ve:
        print(f"Error en los datos de entrada: {ve}")
        return None
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")
        return None

def download_chirps_data(url, output_path):
    """
    Descarga un archivo de datos CHIRPS desde una URL y lo guarda en el disco.
    """
    try:
        response = requests.get(url)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"Datos descargados exitosamente y guardados en {output_path}")
        else:
            print(f"Error al descargar los datos. Código de estado: {response.status_code}")
            return None
        
        try:
            nc_data = Dataset(output_path, mode='r')
            return nc_data
        except Exception as e:
            print(f"Error al cargar el archivo NetCDF: {e}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"Error durante la solicitud de descarga: {e}")
        return None

def download_data_chirps(params):
    base_path = params['save_path']
    chirps_path = os.path.join(base_path, 'chirps')
    os.makedirs(chirps_path, exist_ok=True)
    
    url = generate_chirps_url(params['start_day'], params['start_month'], params['start_year'], params['end_day'], params['end_month'], params['end_year'], params['lat_min'], params['lat_max'], params['lon_min'], params['lon_max'])
    download_chirps_data(url, os.path.join(chirps_path, 'chirps_daily.nc'))
    
    return "Datos CHIRPS descargados exitosamente."
def kendall_tau(a, b):
    """
    Calcula el coeficiente de correlación de Kendall Tau entre dos arrays.

    Parámetros:
    - a (array-like): Primer conjunto de datos.
    - b (array-like): Segundo conjunto de datos.

    Retorna:
    - float: El coeficiente de correlación de Kendall Tau.
    
    Lanza:
    - ValueError: Si las entradas no tienen la misma longitud.
    - Exception: Si ocurre un error inesperado durante el cálculo.
    """
    
    try:
        # Verificar si las entradas tienen la misma longitud
        if len(a) != len(b):
            raise ValueError("Los arrays de entrada deben tener la misma longitud.")
        
        # Calcular el coeficiente de correlación de Kendall Tau
        tau, _ = kendalltau(a, b)
        
        return tau
    
    except ValueError as ve:
        print(f"Error en los datos de entrada: {ve}")
        return None
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")
        return None
def download_data_gcm(params):
    base_path = params['save_path']
    chirps_path = os.path.join(base_path, 'chirps')
    descarga_path = os.path.join(base_path, 'descarga')
    ensamble_path = os.path.join(base_path, 'ensamble')
    os.makedirs(chirps_path, exist_ok=True)
    os.makedirs(descarga_path, exist_ok=True)
    os.makedirs(ensamble_path, exist_ok=True)

    lat_1 = convert_lat_lon(params['lat_1'])
    lat_2 = convert_lat_lon(params['lat_2'])
    lon_1 = convert_lat_lon(params['lon_1'])
    lon_2 = convert_lat_lon(params['lon_2'])
    models = params['models']

    if isinstance(models, str):
        models = [models]

    try:
        urls = list(map(lambda model: get_gcm_url(params['start_time'], params['start_year'], params['last_year'], params['start_step'], params['last_step'], lat_1, lat_2, lon_1, lon_2, model), models))
    except ValueError as e:
        return f"Error al generar URLs: {e}"

    generate_path = lambda model: f"{descarga_path}/{model}_{params['start_time']}_{params['start_step']}-{params['last_step']}.nc"
    paths = list(map(generate_path, models))
    try:
        list(map(download_nc_file, urls, paths))
    except Exception as e:
        return f"Error al descargar archivos NC: {e}"

    start_time = convert_month(params['start_time'])
    start_month = int(start_time + params['start_step'] - 0.5)
    end_month = int(start_time + params['last_step'] - 0.5)
    if end_month > 12:
        end_month -= 12
        years = [params['start_year'], params['start_year'] + 1]
    else:
        years = [params['start_year']]
    months = [start_month, end_month]

    precip_monthly = calculate_winter_precip(os.path.join(chirps_path, 'chirps_daily.nc'), months, params['start_year'], params['last_year'])
    if precip_monthly is None:
        return "Error en el cálculo de precipitación"

    precip_monthly = precip_monthly.rename({'year': 'S'})
    gcm = []
    for file_path in paths:
        try:
            ds = xr.open_dataset(file_path, decode_times=False).aprod
            ds = ds.rename({'L': 'M'})
            units = ds['S'].attrs['units']
            calendar = ds['S'].attrs['calendar']
            origin = pd.Timestamp(units.split('since')[-1].strip())
            dates = [origin + pd.DateOffset(months=month) for month in ds["S"]]
            date_plus = [date + pd.DateOffset(months=math.floor(params['last_step'])) for date in dates]
            years = [date.year for date in date_plus]
            ds_last = ds.assign_coords({'S': years}).copy()
            ds_regrid = xc.regrid(ds_last, precip_monthly.X, precip_monthly.Y)
            gcm.append(ds_regrid)
        except Exception as e:
            return f"Error al procesar archivo {file_path}: {e}"

    try:
        Model_hindcast = xr.concat(gcm, 'M').assign_coords({'M': models})
        Ensamble = Model_hindcast.sel(M=models).mean(dim='M').rename('prec')
        Ensamble.to_netcdf(os.path.join(ensamble_path, 'ensamble.nc'))
    except Exception as e:
        return f"Error al crear el ensamble: {e}"

    return "Datos GCM descargados y procesados exitosamente."

def download_data_gcm_forecast(params):
    base_path = params['save_path']
    ensamble_path = os.path.join(base_path, 'ensamble')
    os.makedirs(ensamble_path, exist_ok=True)

    lat_1 = convert_lat_lon(params['lat_1'])
    lat_2 = convert_lat_lon(params['lat_2'])
    lon_1 = convert_lat_lon(params['lon_1'])
    lon_2 = convert_lat_lon(params['lon_2'])
    models = params['models']

    urls = list(map(lambda model: get_gcm_url(params['start_time'], params['start_year'], params['last_year'], params['start_step'], params['last_step'], lat_1, lat_2, lon_1, lon_2, model), models))
    generate_path = lambda model: f"{ensamble_path}/{model}_{params['start_time']}_{params['start_step']}-{params['last_step']}.nc"
    paths = list(map(generate_path, models))

    list(map(download_nc_file, urls, paths))
    return paths

def start_download_chirps():
    try:
        params = {
            'start_day': entry_start_day_chirps.get(),
            'start_month': entry_start_month_chirps.get(),
            'start_year': int(entry_start_year_chirps.get()),
            'end_day': entry_end_day_chirps.get(),
            'end_month': entry_end_month_chirps.get(),
            'end_year': int(entry_end_year_chirps.get()),
            'lat_min': entry_lat_min_chirps.get(),
            'lat_max': entry_lat_max_chirps.get(),
            'lon_min': entry_lon_min_chirps.get(),
            'lon_max': entry_lon_max_chirps.get(),
            'save_path': entry_save_path_chirps.get()
        }
        result = download_data_chirps(params)
        messagebox.showinfo("Resultado", result)
    except ValueError as ve:
        messagebox.showerror("Error de valor", f"Error en los valores de entrada: {ve}")
    except Exception as e:
        messagebox.showerror("Error", f"Ocurrió un error: {e}")

def start_download_gcm():
    try:
        params = {
            'start_time': entry_start_time_gcm.get(),
            'start_year': int(entry_start_year_gcm.get()),
            'last_year': int(entry_last_year_gcm.get()),
            'start_step': float(entry_start_step_gcm.get()),
            'last_step': float(entry_last_step_gcm.get()),
            'lat_1': entry_lat_1_gcm.get(),
            'lat_2': entry_lat_2_gcm.get(),
            'lon_1': entry_lon_1_gcm.get(),
            'lon_2': entry_lon_2_gcm.get(),
            'models': [models_listbox_gcm.get(idx) for idx in models_listbox_gcm.curselection()],
            'save_path': entry_save_path_chirps.get()  # Reutilizar save_path de CHIRPS
        }
        result = download_data_gcm(params)
        messagebox.showinfo("Resultado", result)
    except ValueError as ve:
        messagebox.showerror("Error de valor", f"Error en los valores de entrada: {ve}")
    except Exception as e:
        messagebox.showerror("Error", f"Ocurrió un error: {e}")

def select_directory(entry_widget):
    path = filedialog.askdirectory()
    if path:
        entry_widget.delete(0, tk.END)
        entry_widget.insert(0, path)

def pearson_to_kendall(r):
    return (2 / np.pi) * np.arcsin(r)

def create_values(limit_A, limit_B, limit_C):
    A = list(range(1, limit_A + 1))
    B = list(range(1, limit_B + 1))
    C = list(range(1, limit_C + 1))
    return A, B, C

def run_forecast_analysis():
    try:
        params = {
            'start_time': entry_start_time_forecast.get(),
            'start_year': int(entry_start_year_gcm.get()),  # Reutilizar start_year de GCM
            'last_year': int(entry_last_year_forecast.get()),
            'start_step': float(entry_start_step_forecast.get()),
            'last_step': float(entry_last_step_forecast.get()),
            'lat_1': entry_lat_1_gcm.get(),  # Reutilizar lat_1 de GCM
            'lat_2': entry_lat_2_gcm.get(),  # Reutilizar lat_2 de GCM
            'lon_1': entry_lon_1_gcm.get(),  # Reutilizar lon_1 de GCM
            'lon_2': entry_lon_2_gcm.get(),  # Reutilizar lon_2 de GCM
            'models': [models_listbox_forecast.get(idx) for idx in models_listbox_forecast.curselection()],
            'save_path': entry_save_path_chirps.get(),  # Reutilizar save_path de CHIRPS
            'eof_x': int(entry_eof_x_forecast.get()),
            'eof_y': int(entry_eof_y_forecast.get()),
            'cca_modos': int(entry_cca_modos_forecast.get())
        }

        base_path = params['save_path']
        chirps_path = os.path.join(base_path, 'chirps')
        ensamble_path = os.path.join(base_path, 'ensamble')
        fore_path = os.path.join(base_path, 'forecast')
        os.makedirs(chirps_path, exist_ok=True)
        os.makedirs(ensamble_path, exist_ok=True)
        os.makedirs(fore_path, exist_ok=True)

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        honduras = world[world.name == 'Honduras']
        current_year = datetime.now().year
        urls = list(map(lambda model: get_gcm_url(params['start_time'], params['start_year'] + (current_year - params['start_year']) - 5, params['last_year'], params['start_step'], params['last_step'], params['lat_1'], params['lat_2'], params['lon_1'], params['lon_2'], model), params['models']))
        generate_path = lambda model: f"{ensamble_path}/{model}_{params['start_time']}_{params['start_step']}-{params['last_step']}.nc"
        paths = list(map(generate_path, params['models']))
        list(map(download_nc_file, urls, paths))

        start_time = convert_month(params['start_time'])
        start_step = params['start_step']
        last_step = params['last_step']
        start_month = int(start_time + start_step - 0.5)
        end_month = int(start_time + last_step - 0.5)
        months = [start_month, end_month]

        obs = calculate_winter_precip(os.path.join(chirps_path, 'chirps_daily.nc'), months, params['start_year'], params['last_year'])
        obs = obs.rename({'year': 'S'})
        drymask = xc.drymask(obs, dry_threshold=0.1, quantile_threshold=0.1)
        obs = obs * drymask

        if obs.rio.crs is None:
            obs = obs.rio.write_crs("EPSG:4326")

        model1_h = xr.open_dataset(os.path.join(ensamble_path, 'ensamble.nc'), decode_times=False).prec
        model1_h = model1_h.expand_dims('M').assign_coords(M=[0])

        if model1_h.rio.crs is None:
            model1_h = model1_h.rio.write_crs("EPSG:4326")

        drymask = xc.drymask(model1_h, dry_threshold=0.1, quantile_threshold=0.1)
        model1_h = model1_h * drymask
        model1_f = model1_h.dropna(dim='S', how='all')
        model1_f_S = model1_f['S'].values
        obs_S = obs['S'].values
        common_S = np.intersect1d(model1_f_S, obs_S)
        model1_f = model1_f.sel(S=common_S)
        obs_f = obs.sel(S=common_S)
        gcm = []
        for file_path in paths:
            try:
                ds = xr.open_dataset(file_path, decode_times=False).aprod
                ds = ds.rename({'L': 'M'})
                units = ds['S'].attrs['units']
                calendar = ds['S'].attrs['calendar']
                origin = pd.Timestamp(units.split('since')[-1].strip())
                dates = [origin + pd.DateOffset(months=month) for month in ds["S"].values]
                date_plus = [date + pd.DateOffset(months=math.floor(last_step)) for date in dates]
                years = [date.year for date in date_plus]
                ds_last = ds.assign_coords({'S': years}).copy()
                ds_regrid = xc.regrid(ds_last, obs_f.X, obs_f.Y)
                gcm.append(ds_regrid)
            except Exception as e:
                print(f"Error al procesar el archivo {file_path}: {e}")

        if gcm:
            Model_hindcast = xr.concat(gcm, 'M')
            print("Concatenación completada exitosamente.")
        else:
            print("No se pudo cargar ningún archivo correctamente.")
            return
        Model_hindcast = Model_hindcast.mean(dim='M')
        fore = Model_hindcast.isel(S=-1, drop=False)

        eof_x = params['eof_x']
        eof_y = params['eof_y']
        cca_modos = params['cca_modos']
        A, B, C = create_values(eof_x, eof_y, cca_modos)
        combinaciones = itertools.product(A, B, C)
        combinaciones_filtradas = [(a, b, c) for a, b, c in combinaciones if c <= min(a, b)]
        resultados = []
        step_difference = params['last_step'] - params['start_step'] + 1
        window = 3 if step_difference == 3 else 2
        S_value = 1 if step_difference == 3 else 0
        for combinacion in combinaciones_filtradas:
            est_kwargs = {'search_override': combinacion, 'latitude_weighting': True}
            ytest_list, preds_list = [], []
            cross_validator = xc.CrossValidator(model1_f, obs_f, window=window)
            for xtrain, ytrain, xtest, ytest in cross_validator:
                reg = xc.CCA(**est_kwargs)
                reg.fit(xtrain, ytrain)
                preds = reg.predict(xtest)
                probs = reg.predict_proba(xtest)
                preds_list.append(preds.isel(S=S_value))
                ytest_list.append(ytest.isel(S=S_value))

            hindcasts_det = xr.concat(preds_list, 'S')
            ytest_concat = xr.concat(ytest_list, 'S')
            pearson = xc.Pearson(hindcasts_det, ytest_concat)
            kendall_dataarray = xr.apply_ufunc(
                pearson_to_kendall,
                pearson,
                vectorize=True,
                dask="parallelized",
                output_dtypes=[float]
            )
            tau_kendall = kendall_dataarray.mean().item()
            resultados.append({'Combinacion': combinacion, 'Goodness Index': tau_kendall})

        resultados_df = pd.DataFrame(resultados).round(2)
        resultados_df.to_excel(os.path.join(fore_path, 'cross_validation.xlsx'), index=False)

        max_tau = resultados_df['Goodness Index'].max()
        max_tau_df = resultados_df[resultados_df['Goodness Index'] == max_tau]
        max_tau_df['complexity'] = max_tau_df['Combinacion'].apply(lambda x: sum(x))
        least_complex_model = max_tau_df.loc[max_tau_df['complexity'].idxmin()]
        best_combination = least_complex_model['Combinacion']

        best_kwargs = {
            'search_override': best_combination,
            'latitude_weighting': True,
        }

        hindcasts_det, hindcasts_prob = [], []
        step_difference = params['last_step'] - params['start_step'] + 1
        window = 3 if step_difference == 3 else 2
        S_value = 1 if step_difference == 3 else 0
        for xtrain, ytrain, xtest, ytest in xc.CrossValidator(model1_f, obs_f, window=window):
            reg = xc.CCA(**best_kwargs)
            reg.fit(xtrain, ytrain)
            preds = reg.predict(xtest)
            probs = reg.predict_proba(xtest)
            hindcasts_det.append(preds.isel(S=S_value))
            hindcasts_prob.append(probs.isel(S=S_value))

        hindcasts_det = xr.concat(hindcasts_det, 'S')
        hindcasts_prob = xr.concat(hindcasts_prob, 'S')

        hindcasts_det = xc.gaussian_smooth(hindcasts_det, kernel=9)
        hindcasts_prob = xc.gaussian_smooth(hindcasts_prob, kernel=9)
        obs = xc.gaussian_smooth(obs, kernel=9)

        if 'M' not in fore.dims:
            fore = fore.expand_dims('M').assign_coords(M=[0])
        if 'S' not in fore.dims:
            fore = fore.expand_dims('S')
        fore['M'] = model1_f['M']
        reg.fit(model1_f, obs_f)
        preds = reg.predict(fore)
        probs = reg.predict_proba(fore)

        forecasts_det_smooth = xc.gaussian_smooth(preds, kernel=9)
        forecasts_prob_smooth = xc.gaussian_smooth(probs, kernel=9)
        forecasts_det_smooth_anomaly = forecasts_det_smooth - hindcasts_det.mean('S')

        forecasts_det_smooth.to_netcdf(os.path.join(fore_path, 'forecasts_deterministic.nc'))
        forecasts_det_smooth_anomaly.to_netcdf(os.path.join(fore_path, 'forecasts_anomaly.nc'))
        forecasts_prob_smooth.to_netcdf(os.path.join(fore_path, 'forecasts_tercile.nc'))
        forecsts_probs20 = reg.predict_proba(fore, quantile=0.2)
        forecsts_probs80 = 1 - reg.predict_proba(fore, quantile=0.8)
        forecasts_prob20_smooth = xc.gaussian_smooth(forecsts_probs20, kernel=9)
        forecasts_prob80_smooth = xc.gaussian_smooth(forecsts_probs80, kernel=9)
        forecasts_prob20_smooth.to_netcdf(os.path.join(fore_path, 'forecasts_prob_20.nc'))
        forecasts_prob80_smooth.to_netcdf(os.path.join(fore_path, 'forecasts_prob_80.nc'))

        # Mostrar el DataFrame con los valores de validación
        print(resultados_df)
        
        # Mostrar el DataFrame en una ventana emergente
        display_dataframe(resultados_df)

        # Mostrar un mensaje con la mejor combinación
        messagebox.showinfo("Mejor Combinación", f"La mejor combinación es: {best_combination}")

    except ValueError as ve:
        messagebox.showerror("Error de valor", f"Error en los valores de entrada: {ve}")
    except Exception as e:
        messagebox.showerror("Error", f"Ocurrió un error: {e}")

def display_dataframe(df):
    """
    Muestra el DataFrame en una ventana emergente.
    """
    top = tk.Toplevel()
    top.title("Resultados de Validación")
    frame = ttk.Frame(top)
    frame.pack(fill='both', expand=True)
    pt = ttk.Treeview(frame)
    pt.pack(fill='both', expand=True)

    pt["columns"] = list(df.columns)
    pt["show"] = "headings"
    for col in df.columns:
        pt.heading(col, text=col)
        pt.column(col, width=100)
    for index, row in df.iterrows():
        pt.insert("", "end", values=list(row))

    scrollbar = ttk.Scrollbar(frame, orient="vertical", command=pt.yview)
    pt.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side='right', fill='y')

# Crear la ventana principal
root = tk.Tk()
root.title("Descarga y Procesamiento de Datos de Precipitación")
root.geometry("600x800")
root.configure(bg='white')

# Estilo de la aplicación
style = ttk.Style()
style.configure("TLabel", font=("Arial", 12), background="white")
style.configure("TButton", font=("Arial", 12), padding=6)
style.configure("TEntry", font=("Arial", 12), padding=6)
style.configure("TNotebook.Tab", font=("Arial", 12), padding=6)

# Crear el cuaderno de pestañas
notebook = ttk.Notebook(root)
notebook.pack(pady=10, expand=True)

# Crear los marcos para cada pestaña
frame_chirps = ttk.Frame(notebook, style="TFrame", padding=10)
frame_gcm = ttk.Frame(notebook, style="TFrame", padding=10)
frame_forecast = ttk.Frame(notebook, style="TFrame", padding=10)

# Agregar las pestañas al cuaderno
notebook.add(frame_chirps, text='Datos CHIRPS')
notebook.add(frame_gcm, text='Modelos GCM')
notebook.add(frame_forecast, text='Forecast')

def create_label_entry(parent, label_text, default_value):
    frame = ttk.Frame(parent)
    frame.pack(pady=5, padx=10, fill='x')
    label = ttk.Label(frame, text=label_text, width=20)
    label.pack(side='left')
    entry = ttk.Entry(frame)
    entry.pack(side='left', fill='x', expand=True)
    entry.insert(0, default_value)
    return entry

# Contenido de la pestaña CHIRPS
ttk.Label(frame_chirps, text="Parámetros CHIRPS", font=("Arial", 14, "bold")).pack(pady=10)
entry_start_day_chirps = create_label_entry(frame_chirps, "Start Day", "01")
entry_start_month_chirps = create_label_entry(frame_chirps, "Start Month", "Jan")
entry_start_year_chirps = create_label_entry(frame_chirps, "Start Year", "1982")
entry_end_day_chirps = create_label_entry(frame_chirps, "End Day", "31")
entry_end_month_chirps = create_label_entry(frame_chirps, "End Month", "Dec")
entry_end_year_chirps = create_label_entry(frame_chirps, "End Year", "2023")
entry_lat_min_chirps = create_label_entry(frame_chirps, "Latitud Mínima", "12N")
entry_lat_max_chirps = create_label_entry(frame_chirps, "Latitud Máxima", "16N")
entry_lon_min_chirps = create_label_entry(frame_chirps, "Longitud Mínima", "90W")
entry_lon_max_chirps = create_label_entry(frame_chirps, "Longitud Máxima", "83W")
entry_save_path_chirps = create_label_entry(frame_chirps, "Ruta de Guardado", "C:/ruta/de/descarga")
ttk.Button(frame_chirps, text="Seleccionar Directorio", command=lambda: select_directory(entry_save_path_chirps)).pack(pady=5)
ttk.Button(frame_chirps, text="Descargar Datos CHIRPS", command=start_download_chirps).pack(pady=5)
status_chirps = tk.StringVar()
ttk.Label(frame_chirps, textvariable=status_chirps, background="white").pack(pady=5)

# Contenido de la pestaña GCM
ttk.Label(frame_gcm, text="Parámetros GCM", font=("Arial", 14, "bold")).pack(pady=10)
entry_start_time_gcm = create_label_entry(frame_gcm, "Start Time", "Mar")
entry_start_year_gcm = create_label_entry(frame_gcm, "Start Year", "1982")
entry_last_year_gcm = create_label_entry(frame_gcm, "Last Year", "2023")
entry_start_step_gcm = create_label_entry(frame_gcm, "Start Step", "5.5")
entry_last_step_gcm = create_label_entry(frame_gcm, "Last Step", "7.5")
entry_lat_1_gcm = create_label_entry(frame_gcm, "Latitud 1", "12N")
entry_lat_2_gcm = create_label_entry(frame_gcm, "Latitud 2", "16N")
entry_lon_1_gcm = create_label_entry(frame_gcm, "Longitud 1", "90W")
entry_lon_2_gcm = create_label_entry(frame_gcm, "Longitud 2", "83W")
ttk.Label(frame_gcm, text="Modelos").pack(pady=5)
models_listbox_gcm = tk.Listbox(frame_gcm, selectmode=tk.MULTIPLE, font=("Arial", 12), height=6)
models = ["NCEP-CFSv2", "Can-1", "Can-2", "CCSM4", "GFDL", "NASA", 
          "ECMWF", "Meteo-France", "Glosea", "CMCC", "DWD", "JMA"]
for model in models:
    models_listbox_gcm.insert(tk.END, model)
models_listbox_gcm.pack(pady=5)
ttk.Button(frame_gcm, text="Descargar Datos GCM", command=start_download_gcm).pack(pady=5)
status_gcm = tk.StringVar()
ttk.Label(frame_gcm, textvariable=status_gcm, background="white").pack(pady=5)

# Contenido de la pestaña Forecast
ttk.Label(frame_forecast, text="Parámetros Forecast", font=("Arial", 14, "bold")).pack(pady=10)
canvas_forecast = tk.Canvas(frame_forecast, bg='white')
scroll_y_forecast = ttk.Scrollbar(frame_forecast, orient="vertical", command=canvas_forecast.yview)
frame_forecast_content = ttk.Frame(canvas_forecast, style="TFrame", padding=10)
canvas_forecast.create_window((0, 0), window=frame_forecast_content, anchor='nw')
canvas_forecast.update_idletasks()
canvas_forecast.configure(scrollregion=canvas_forecast.bbox('all'), yscrollcommand=scroll_y_forecast.set)
canvas_forecast.pack(side='left', fill='both', expand=True, padx=10, pady=10)
scroll_y_forecast.pack(side='right', fill='y')
entry_start_time_forecast = create_label_entry(frame_forecast_content, "Start Time", "Mar")
entry_last_year_forecast = create_label_entry(frame_forecast_content, "Forecast Year", "2024")
entry_start_step_forecast = create_label_entry(frame_forecast_content, "Start Step", "5.5")
entry_last_step_forecast = create_label_entry(frame_forecast_content, "Last Step", "7.5")
ttk.Label(frame_forecast_content, text="Modelos").pack(pady=5)
models_listbox_forecast = tk.Listbox(frame_forecast_content, selectmode=tk.MULTIPLE, font=("Arial", 12), height=6)
for model in models:
    models_listbox_forecast.insert(tk.END, model)
models_listbox_forecast.pack(pady=5)
entry_eof_x_forecast = create_label_entry(frame_forecast_content, "EOF X", "1")
entry_eof_y_forecast = create_label_entry(frame_forecast_content, "EOF Y", "1")
entry_cca_modos_forecast = create_label_entry(frame_forecast_content, "CCA Modos", "1")
ttk.Button(frame_forecast_content, text="Descargar y Analizar Forecast", command=run_forecast_analysis).pack(pady=5)
status_forecast = tk.StringVar()
ttk.Label(frame_forecast_content, textvariable=status_forecast, background="white").pack(pady=5)

frame_forecast_content.bind("<Configure>", lambda e: canvas_forecast.configure(scrollregion=canvas_forecast.bbox("all")))

# Iniciar la aplicación
root.mainloop()