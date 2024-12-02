"""
Este programa identifica las secuencias de exones de un gen en el genoma de cloroplasto de
Solanum phureja. Utiliza dos archivos en formato FASTA: uno contiene un único contig (Genoma)
y el otro las secuencias de los exones, ambos en formato de nucleótidos. En los exones las
líneas de encabezado deben de tener el indicador del exon deben de estar en la tercera posicion,
se recomienda comenzar con >, seguidas del nombre del gen (sin espacios) y  identificador del
exón (también sin espacios).Ejm: > PDS exon1. Ademas estos deben de ser ordenados segun su lugar
en el genoma
"""

# Importar herramientas
import os

# Funciones
def complement(genome):
    complementary= {"A":"T", "C":"G", "G":"C", "T":"A"}
    c_genome =""
    for letter in genome:
        c_genome += complementary[letter]
    return c_genome

def reverse(genome):
    return genome[::-1]

def reverse_complement(genome):
    c_genome = complement(genome)
    reverse_complement_genome= reverse(c_genome)
    return reverse_complement_genome

# Variables iniciales
genome = str()
reversecomplem_genome = str()
exons = {}
num_exon = str()
seq = str()
positions = {}
exons_starts1 = []
exons_starts2 = []
list_positions = []
min_difference = float("inf")
position_starts_exon1 = int()
position_starts_exon2 = int()
gene_name = str()
genome_name = str()
max_intron_length = 3000

# Primero ingreso el genoma en donde voy a buscar el gen
with open("Solanum_tuberosum_phureja_genome_plastid_nucl.fasta") as file:
    for line in file:
# Proceso el genoma para que sea mas facil trabajar con el, dejandolo en la misma linea
        if line.startswith(">"):
            pass
        else:
            line = line.rstrip("\n")
            genome += line
genome = genome.upper() # Asegurarme que todos los nucleotidos esten en mayusculas

# Creo el reverso complementerio del genoma
reversecomplem_genome = reverse_complement(genome)

# Ingreso el archivo con los exones que voy a buscar
with open("exon_petB_nucl.fasta") as f:
    for l in f:
        if l.startswith(">"):
            if exons:
                exons[num_exon] = seq
            temp = l.split()  # Variable temporal para el titulo
            num_exon = temp[2]
            gene_name = temp[1]
            seq = "" # Poner en blanco la variable para que la siguiente corrida no tenga la secuencia del exon anterior
        else:
            seq += l.strip()
            exons[num_exon] = seq

# Ahora voy a buscar todas las posiciones donde estan los exones dentro del genoma
try: # Correr en direccion forward
    for exon_id, exon_seq in exons.items():
        start = 0
        list_positions = [] # Resetear la lista para cada exon
        while True:
            start = genome.find(exon_seq, start)
            if start == -1:
                break
            start += 1
            end = start + len(exon_seq) - 1
            list_positions.append({"start": start, "end": end})
            start += len(exon_seq)
            positions[exon_id] = list_positions # Cada exon tiene todas sus posibles posiciones
            seq_genome = genome
    if not positions:
        raise ValueError()
except: # Correr en direccion reverse
    for exon_id, exon_seq in exons.items():
        start = 0
        list_positions = [] # Resetear la lista para cada exon
        while True:
            start = reversecomplem_genome.find(exon_seq, start)
            if start == -1:
                break
            start += 1
            end = start + len(exon_seq) - 1
            list_positions.append({"start": start, "end": end})
            start += len(exon_seq)
            positions[exon_id] = list_positions
            seq_genome = reversecomplem_genome

for exon_id, positions_exon in positions.items():
    for pos in positions_exon:
        if exon_id == "exon1":
            exons_starts1.append(pos["start"]) # Pongo las posiciones del exon 1 en solo una lista
        elif exon_id == "exon2":
            exons_starts2.append(pos["start"]) # Pongo las posiciones del exon 2 en solo una lista

for e in range(len(exons_starts1)):
    dif= exons_starts2[0] - exons_starts1[e] # Porque el exon 2 solo tiene una posicion posible
    if 0 < dif < max_intron_length and dif < min_difference:
        min_difference = dif
        position_starts_exon1 = exons_starts1[e]
        position_starts_exon2 = exons_starts2[0]

# Organizo los datos para el archivo de salida
genome_length = len(genome)
if seq_genome == reversecomplem_genome:
    forward_start_exon1 = genome_length - position_starts_exon2 - len(exons["exon2"]) + 2
    forward_end_exon1 = genome_length - position_starts_exon2 + 1
    exon1_seq = genome[forward_start_exon1 : forward_end_exon1]
    forward_start_exon2 = genome_length - position_starts_exon1 - len(exons["exon1"])+ 2
    forward_end_exon2 = genome_length - position_starts_exon1 + 1
    exon2_seq = genome[forward_start_exon2: forward_end_exon2]
    intron1_seq = genome[forward_end_exon1 : forward_start_exon2]
    gene_seq = genome[forward_end_exon2 : forward_start_exon1]
    position_starts_exon1 = forward_start_exon1
    position_starts_exon2 = forward_start_exon2
    position_end_exon1 = forward_end_exon1
    position_end_exon2 = forward_end_exon2
else:
    # Calcular las posiciones directamente
    position_end_exon1 = position_starts_exon1 + len(exons["exon1"]) - 1
    position_end_exon2 = position_starts_exon2 + len(exons["exon2"]) - 1
    exon1_seq = seq_genome[position_starts_exon1 - 1: position_end_exon1]
    intron1_seq = seq_genome[position_end_exon1: position_starts_exon2 - 1]
    exon2_seq = seq_genome[position_starts_exon2 - 1: position_end_exon2]
    gene_seq = seq_genome[position_starts_exon1 - 1: position_end_exon2]

# Calcular el reverso complementario del gen
rc_gene_seq = reverse_complement(gene_seq)

# Generar el documento de salida
genome_name = "Solanum_tuberosum_phureja_genome_plastid_nucl.fasta"
genome_name = genome_name.split("_")
genome_name = genome_name[0] + " "+genome_name[1] + " " + genome_name[2]
source = "new_file.txt"  # Nombre del archivo temporal

with open(source, 'w') as fl:
    fl.write(f"El gen {gene_name} en el genoma de {genome_name} tiene las sigueintes caracateristicas:\n"
             f"\nEl primer exon va de la posicion {position_starts_exon1} a la posicion {position_end_exon1} dentro del genoma. Con una secuencia: \n{exon1_seq}\n"
             f"\nEl primer intron va de la posicion {position_end_exon1} a la posicion {position_starts_exon2} dentro del genoma. Con una secuencia: \n{intron1_seq}\n"
             f"\nEl segundo exon va de la posicion {position_starts_exon2} a la posicion {position_end_exon2} dentro del genoma. Con una secuencia: \n{exon2_seq}\n"
             f"\nAsi que el gen completo va de la posicion {position_starts_exon1} a la posicion {position_end_exon2} dentro del genoma. Con una secuencia: \n{gene_seq}\n"
             f"\nY la secuenci del reverso complementario de este gen es: \n{rc_gene_seq}")
dest = f"{gene_name}_exons_positions_genome.txt" # Nuevo nombre del documento de salida
os.rename(source, dest) # Renombro el documento de salida