#!/usr/bin/python3

# Este Script encuentra anticodones flanqueados por secuencias palindromas.
# El resultado es la impresion en pantalla de un archivo fasta con seceuncias candidatas a ser tRNAs

from Bio import SeqIO
import re
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: tRNAdetector.py <archivo.fasta> <Codon>")
        sys.exit(1)

archivo_fasta =open(sys.argv[1])
codon =(sys.argv[2])

def reverse_complement(secuencia):
    # Calcula la reversa complementaria de una secuencia

    complementarias = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversa_complementaria = ''.join(complementarias.get(base, base) for base in secuencia[::-1])  # Reemplazo manual y reverso
    return reversa_complementaria


def es_palindromo_complementario(secuencia):
  #Verifica si una secuencia de ADN es palindrómica complementaria.

  secuencia_complementaria = reverse_complement(secuencia)
  return secuencia == secuencia_complementaria



# Cargamos las secuencias del archivo FASTA
secuencias = list(SeqIO.parse(archivo_fasta, "fasta"))

def palindromo_acceptor_setm(s5,s3):
    seq=(s5,s3)
    print(seq)

def buscar_anticodon_loop(secuencias,motivo):
    #Inicializa variable para ids del fasta
    id = 1
    for secuencias in secuencias:
        coincidencias = re.finditer(str(motivo), str(secuencias.seq))
        if coincidencias:
            for coincidencia in coincidencias:
                inicio = coincidencia.start()
                fin = coincidencia.end()
                # 5 indica la distancia al final del palindromo despues del anticodon; 2 indica que hay 2 nucleootidos despues del anticodon que no son tomados en cuenta (pero forman parte del loop del brazo del anticodon).
                # en este ejemplo el palindromo resultante es de 3 nucleotidos
                anticodon_loop = secuencias.seq[inicio - 5:inicio-2] + secuencias.seq[fin+2:fin + 5]
                acceptor_setm5prime = secuencias.seq[inicio-30:inicio-20]
                acceptor_setm3prime = secuencias.seq[fin+26:fin+30]
                if es_palindromo_complementario(anticodon_loop):
                    # En un futuro implementar la detección del D-arm y TψC-arm
                    # for i in range(1,len(acceptor_setm5prime)):
                    #     s5 = (acceptor_setm5prime[0:i+1])
                    #     for j in range(1,len(acceptor_setm3prime)):
                    #         s3 = (acceptor_setm3prime[0:j+1])
                    #         palindromo_acceptor_setm(s5,s3)
                            #imprime la secuencia candidata a ser un tRNA
                    print('>',motivo,id, sep="")
                    print(secuencias.seq[inicio-30:fin+30])
                    id = id +1
anticodon = str(reverse_complement(codon))
motivos = [codon,anticodon]
for motivo in motivos:
    buscar_anticodon_loop(secuencias,motivo)
