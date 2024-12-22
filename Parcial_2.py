import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF
from Bio.Seq import Seq
from Bio import SeqIO
from textwrap import wrap

# Función para ajustar texto en el PDF
def wrap_text(pdf, text, max_width):
    """
    Ajusta el texto al ancho máximo especificado para evitar que se desborde en el PDF.
    """
    lines = wrap(text, width=max_width // 3)  # Ajusta el factor según el tamaño de letra
    for line in lines:
        pdf.cell(0, 10, line, ln=True)

# Función base_frequencies para calcular frecuencias de bases
def base_frequencies(sequence):
    """
    Calcula la frecuencia de las bases nitrogenadas (A, T, G, C) en una secuencia de ADN.
    """
    sequence = sequence.upper()
    bases = ['A', 'T', 'G', 'C']
    frequencies = {base: sequence.count(base) for base in bases}
    return frequencies

# Nueva función para graficar frecuencias de bases de manera más compleja
def plot_base_frequencies(sequence):
    """
    Calcula y grafica la frecuencia de las bases nitrogenadas (A, T, G, C).
    """
    frequencies = base_frequencies(sequence)
    df = pd.DataFrame.from_dict(frequencies, orient='index', columns=['Conteo'])
    df['Frecuencia (%)'] = (df['Conteo'] / len(sequence)) * 100
    df['Base'] = df.index

    plt.figure(figsize=(10, 6))
    sns.barplot(x='Base', y='Frecuencia (%)', data=df, palette='viridis', edgecolor='black')
    plt.title("Frecuencia de bases nitrogenadas", fontsize=16, color="darkblue")
    plt.xlabel("Base", fontsize=14)
    plt.ylabel("Frecuencia (%)", fontsize=14)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    for i, row in df.iterrows():
        plt.text(row.name, row['Frecuencia (%)'] + 0.5, f"{row['Frecuencia (%)']:.2f}%", ha='center', fontsize=12)
    plt.tight_layout()
    plt.savefig("frecuencia_bases_compleja.png")
    plt.close()

# Función para detectar regiones de fagos con una gráfica mejorada
def detect_phage_regions(sequence, window_size, gc_threshold):
    """
    Detecta regiones de fagos en la secuencia de ADN en función del contenido de GC.
    """
    regions = []
    gc_content_list = []
    positions = []

    for i in range(0, len(sequence) - window_size + 1, window_size):
        subsequence = sequence[i:i + window_size]
        subseq_gc_content = gc_content(subsequence)
        gc_content_list.append(subseq_gc_content)
        positions.append(i)

        if subseq_gc_content < gc_threshold:  # Regiones con contenido GC bajo
            regions.append((i, i + window_size, subseq_gc_content))

    gc_global = sum(gc_content_list) / len(gc_content_list) if gc_content_list else 0

    # Gráfica mejorada
    plt.figure(figsize=(12, 6))
    plt.plot(positions, gc_content_list, label="Contenido GC (%)", color="blue", lw=2)
    plt.axhline(gc_global, color="green", linestyle="--", label=f"GC Global ({gc_global:.2f}%)")
    plt.axhline(gc_threshold, color="red", linestyle="--", label=f"Umbral GC ({gc_threshold}%)")
    plt.fill_between(positions, gc_content_list, gc_threshold, where=[gc < gc_threshold for gc in gc_content_list],
                     color="red", alpha=0.3, label="Regiones de fagos")
    plt.title("Detección de regiones de fagos basada en contenido GC", fontsize=16, color="darkred")
    plt.xlabel("Posición en la secuencia", fontsize=14)
    plt.ylabel("Contenido GC (%)", fontsize=14)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.tight_layout()
    plt.savefig("detectar_fagos.png")
    plt.close()

    return regions, gc_global

# Función para calcular el contenido de GC
def gc_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

# Función para calcular frecuencias de codones
def codon_frequencies(sequence):
    seq = Seq(sequence[:len(sequence) - len(sequence) % 3])  # Recortar a múltiplo de 3
    codons = [str(seq[i:i + 3]) for i in range(0, len(seq), 3)]
    codon_counts = pd.Series(codons).value_counts()
    return codon_counts

# Función para encontrar ORFs y calcular sus longitudes
def find_orfs_lengths(sequence):
    seq = Seq(sequence)
    orfs = seq.split("ATG")  # Dividir por posibles inicios de ORFs
    lengths = []
    for orf in orfs[1:]:  # Ignorar el primer fragmento (antes del primer ATG)
        stop_codon = orf.find("TAA") or orf.find("TAG") or orf.find("TGA")
        if stop_codon != -1:
            lengths.append((stop_codon + 3) // 3)  # Longitud en codones
    return lengths

# Nueva Función para graficar la distribución de las bases en segmentos de la secuencia
def plot_base_distribution(sequence, window_size=1000):
    """
    Calcula y grafica la distribución de las bases (A, T, G, C) a lo largo de la secuencia en segmentos.
    """
    bases = ['A', 'T', 'G', 'C']
    base_freqs = {base: [] for base in bases}

    # Calcular frecuencias de las bases en ventanas
    for i in range(0, len(sequence) - window_size + 1, window_size):
        subsequence = sequence[i:i + window_size]
        for base in bases:
            base_freqs[base].append(subsequence.count(base) / window_size)

    # Convertir en DataFrame para facilidad de manejo
    base_freq_df = pd.DataFrame(base_freqs)

    # Graficar la distribución de las bases
    plt.figure(figsize=(10, 6))
    base_freq_df.plot(kind='line', lw=2)
    plt.title("Distribución de bases (A, T, G, C) a lo largo de la secuencia", color="darkred")
    plt.xlabel("Posición de la ventana", fontsize=12, fontweight='bold')
    plt.ylabel("Frecuencia de base", fontsize=12, fontweight='bold')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("distribucion_bases.png")
    plt.close()

# Función principal para crear todos los gráficos y generar el PDF
def create_pdf_report(sequence, filename="reporte.pdf"):
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Helvetica", size=16)

    # 1. Frecuencia de bases (gráfica compleja)
    pdf.cell(0, 10, "Frecuencia de bases (A, T, G, C):", ln=True)
    wrap_text(pdf, "Este gráfico muestra la frecuencia relativa de cada base nitrogenada.", 190)
    plot_base_frequencies(sequence)
    pdf.image("frecuencia_bases_compleja.png", x=10, y=60, w=190)
    pdf.add_page()

    # 2. Detección de fagos
    regions, gc_global = detect_phage_regions(sequence, window_size=1000, gc_threshold=10)
    pdf.cell(0, 10, "Detección de posibles fagos:", ln=True)
    wrap_text(pdf,
              f"Este gráfico muestra regiones con contenido GC desviado en más del 10% del promedio ({gc_global:.2f}%).",
              190)
    pdf.image("detectar_fagos.png", x=10, y=60, w=190)
    pdf.add_page()

    # 3. Mapa de calor del contenido GC
    pdf.cell(0, 10, "Mapa de calor del contenido GC:", ln=True)
    wrap_text(pdf, "Este mapa muestra la distribución del contenido GC en la secuencia analizada.", 190)
    window_size = 100
    gc_values = [gc_content(sequence[i:i + window_size]) for i in
                 range(0, len(sequence) - window_size + 1, window_size)]
    heatmap_data = pd.DataFrame({'Posición': range(0, len(gc_values)), 'Contenido GC (%)': gc_values})
    heatmap_data = heatmap_data.pivot_table(index='Posición', values='Contenido GC (%)')
    plt.figure(figsize=(10, 6))
    sns.heatmap(heatmap_data, cmap="coolwarm", cbar_kws={"label": "Contenido GC (%)"})
    plt.title("Mapa de calor del contenido GC", color="purple")
    plt.tight_layout()
    plt.savefig("mapa_gc.png")
    pdf.image("mapa_gc.png", x=10, y=60, w=190)
    pdf.add_page()

    # 4. Frecuencia de codones
    codon_freq = codon_frequencies(sequence).head(10)
    pdf.cell(0, 10, "Frecuencia de los 10 codones más comunes:", ln=True)
    wrap_text(pdf, "Este gráfico muestra los codones más frecuentes en la secuencia analizada.", 190)
    plt.figure(figsize=(10, 6))
    sns.barplot(x=codon_freq.index, y=codon_freq.values, palette="viridis")
    plt.title("Frecuencia de los 10 codones más comunes", color="darkorange")
    plt.xlabel("Codón")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plt.savefig("frecuencia_codones.png")
    pdf.image("frecuencia_codones.png", x=10, y=60, w=190)
    pdf.add_page()

    # 5. Longitudes de ORFs (distribución)
    orf_lengths = find_orfs_lengths(sequence)
    pdf.cell(0, 10, "Distribución de las longitudes de ORFs:", ln=True)
    wrap_text(pdf, "Este gráfico muestra la distribución de las longitudes de los marcos abiertos de lectura (ORFs).",
              100)
    plt.figure(figsize=(10, 6))
    sns.histplot(orf_lengths, kde=True, bins=20, color="purple")
    plt.title("Distribución de longitudes de ORFs", color="brown")
    plt.xlabel("Longitud (codones)")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plt.savefig("longitudes_orfs.png")
    pdf.image("longitudes_orfs.png", x=10, y=60, w=190)
    pdf.add_page()

    # 6. Distribución de bases en segmentos
    plot_base_distribution(sequence)
    pdf.cell(0, 10, "Distribución de bases en la secuencia:", ln=True)
    wrap_text(pdf, "Este gráfico muestra cómo varían las frecuencias de las bases a lo largo de la secuencia.", 220)
    pdf.image("distribucion_bases.png", x=10, y=60, w=190)
    pdf.add_page()

    # Guardar PDF
    pdf.output(filename)
    print(f"Reporte guardado como {filename}")

# Leer archivo FASTA y generar el informe
fasta_file = "E.coli.fasta"
sequence = ""
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence += str(record.seq)

create_pdf_report(sequence)
