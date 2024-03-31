import pandas as pd
import numpy as np
import csv
import numpy as np
from scipy import stats
"""Librerías de química"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolFromInchi
from rdkit.Chem.inchi import MolToInchi
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import matplotlib.pyplot as plt
import numpy as np                #!pip install numpy
import pandas as pd               #!pip install pandas
import matplotlib.pyplot as plt   #!pip install matplotlib
import seaborn as sns
import rdkit
from rdkit.Chem import Descriptors
import seaborn as sns
from rdkit.DataManip.Metric import GetTanimotoDistMat
from rdkit.DataManip.Metric import GetTanimotoSimMat
from rdkit import rdBase
from rdkit.Chem import RDConfig
import matplotlib.pyplot as plt
from kneed import KneeLocator
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
import pubchempy


"SE SE VA A VERIFICAR QUE LAS SUSTANCIAS Y SUS SMILES CORRESPONDEN"
#Las sustancias aquí no están catalogadas por "activas" o "inactivas"

smiles_all = pd.read_csv('C:/Users/licit/Downloads/ENSAYOS TOXCAST PARA TESIS/smiles_all_concatened.csv', sep=';')
print(smiles_all.columns)


"IMAGENES DE TODAS LAS MOLÉCULAS"
from rdkit.Chem import Draw

images_smiles_all = [] #LISTA DE TODAS LAS IMAGENES DE LAS MOLÉCULAS
for i in smiles_all['SMILES']:
    try:
        m = Chem.MolFromSmiles(i)
        img = Draw.MolToImage(m)
        #images_smiles_all.append(img)
    except TypeError:
        img = 'NA'  # Si se produce un error, asignamos 'NA'
    images_smiles_all.append(img)

#print(images_smiles_all)

"""ELEMENTOS QUE GENERAN NA"""
failed_indices = [i for i, img in enumerate(images_smiles_all) if img == 'nan']

print("Índices de elementos en images_smiles_all que generaron 'NA':", failed_indices)


"""NOMBRE DE LAS MOLECULAS"""

moleculas_name = smiles_all['PREFERRED_NAME'].tolist() #lista de nombres de moléculas

"""TRANSFORMAR LOS NOMBRES EN IMÁGENES"""

from PIL import Image, ImageDraw, ImageFont

def texto_a_imagen(texto, ancho, alto, fuente):
    imagen = Image.new('RGB', (ancho, alto), color='white')
    dibujo = ImageDraw.Draw(imagen)
    dibujo.text((0, 0), texto, fill='black', font=fuente)
    return imagen

# Definir una fuente para el texto
fuente = ImageFont.truetype("arial.ttf", 9)  # Puedes ajustar el tamaño y el tipo de fuente según tus preferencias

# Crear imágenes para cada nombre de molécula
imagenes_nombres = []
for nombre in moleculas_name:
    try:
        imagen_nombre = texto_a_imagen(nombre, 100, 20, fuente)  # Ancho y alto ajustables según sea necesario
    except (AttributeError,TypeError):
        imagen_nombre = 'NA'  # Si se produce un error, asignamos 'NA'
    imagenes_nombres.append(imagen_nombre)

"""ELEMENTOS QUE GENERAN NA"""
failed_indices = [i for i, imagen_nombre in enumerate(imagenes_nombres) if imagen_nombre == 'NA']

print("Índices de elementos en moleculas_name que generaron 'NA':", failed_indices)

#print(imagenes_nombres)

"""CREAR EL MOSAICO PARA LAS MOLÉCULAS Y SUS NOMBRES"""

# Definir el número de columnas para el mosaico
num_columnas = 10 # Puedes ajustar este valor según tus preferencias

# Calcular el número de filas necesario para mostrar todas las imágenes
num_imagenes = len(imagenes_nombres)
num_filas = (num_imagenes - 1) // num_columnas + 1

# Calcular el tamaño total del mosaico
alto_maximo = max(
    imagen.height for imagen in images_smiles_all) + 20  # Agregar espacio para el nombre debajo de cada imagen
ancho_total = num_columnas * max(imagen.width for imagen in images_smiles_all)

# Crear una nueva imagen para el mosaico
mosaico = Image.new('RGB', (ancho_total, alto_maximo * num_filas), color='white')

# Zipear las listas de imágenes de moléculas y nombres
imagenes_zip = zip(images_smiles_all, imagenes_nombres)

# Pegar cada imagen y su nombre en el mosaico
for i, (imagen_molecula, imagen_nombre) in enumerate(imagenes_zip):
    # Calcular las coordenadas para pegar la imagen actual
    fila = i // num_columnas
    columna = i % num_columnas
    x_offset = columna * imagen_molecula.width
    y_offset = fila * alto_maximo

    # Pegar la imagen de la molécula
    mosaico.paste(imagen_molecula, (x_offset, y_offset))

    # Pegar la imagen del nombre debajo de la imagen de la molécula
    mosaico.paste(imagen_nombre, (x_offset, y_offset + imagen_molecula.height))

# Mostrar el mosaico
#mosaico.show()
# Guardar el mosaico como una imagen JPEG con una compresión ajustada y una resolución reducida
#mosaico.save("mosaico_moleculas.jpg")  # Adjust quality as needed

"OBTENCIÓN DE LOS IUPAC FROM SMILES"

IUPAC_NAME = []
for i in smiles_all['SMILES']:
    try:
        compounds = pubchempy.get_compounds(i, namespace='smiles')
        match = compounds[0]
        a = match.iupac_name
    except pubchempy.BadRequestError:
        a = 'NA'  # Si se produce un error, asignamos 'NA'
    print(a)
    IUPAC_NAME.append(a)

smiles_all['IPUAC'] = IUPAC_NAME
smiles_all.to_csv('smiles_iupac.csv')



