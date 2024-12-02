# Gene_identification
Este es un programa para identificar exones e intrones de un gen especifico dentro de un genoma

## Archivos de Entrada
El programa espera dos archivos en formato fasta.
Un archivo de entrada que contiene la secuencia de nucleotidos del genoma en el que se va a buscar el gen en un solo conting. Se recomienda que la descripcion del organismo este en el nombre del archivo separado con "_" y que la informacion importante no supere los 3 primeros items, como se muestra en el genoma de prueba "Solanum_tuberosum_phureja_genome_plastid_nucl.fasta"
Un archivo de entrada que contiene las secuencias de nucleotidos de los exones del gen que se desea encontrar cada una identificada por el formato fasta. Se recomienda que el titulo de cada conting tenga el formato > Nombre_del_gen Indicativo_del_exon

## Archivos de Salida
La salida es un archivo de texto, en donde estara la posicion de los intrones y exones con sus respectivas secuencias, ademas de la secuencia de la posicion y secuencia del gen dentro.
