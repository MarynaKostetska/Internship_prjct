from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import logging
import os

""" Налаштування логування
Логування дозволяє отримувати додаткову інформацію про виконання коду,
наприклад, чи було знайдено промотори або чи є проблеми з вхідними файлами."""

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def find_promotors_with_parser(gff_file, fasta_file, result_folder):
    # Визначення назви файлу без розширення
    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]

    # Динамічне створення імен вихідних файлів
    promoter_output = os.path.join(result_folder, f"{fasta_basename}_promoters.fasta")
    non_promoter_output = os.path.join(result_folder, f"{fasta_basename}_non_promoters.fasta")

    # Створення вихідної папки, якщо вона ще не існує
    os.makedirs(result_folder, exist_ok = True)

    #Read FASTA
    record = next(SeqIO.parse(fasta_file, "fasta"))

    #Read GFF
    with open(gff_file) as gff_handle:
        gff_iterator = GFF.parse(gff_handle, base_dict = {record.id: record})
        genome_record = next(gff_iterator) #should be the same ID like in FASTA

    #Collect promoters and non promoters
    promoters = []
    non_promoters = []

    for feature in genome_record.features:
        # Перевірка типу feature
        if feature.type == "regulatory" and any("promoter" in note.lower() for note in feature.qualifiers.get("note", [])):
            promoters.append(feature.location)

        else:
            non_promoters.append(feature.location)

    # Check if no promoters were found
    if not promoters:
        logger.warning("No promoters found in the GFF file.")

    #Save sequences separately
    with open(promoter_output, "w") as prom_file:
        for i, location in enumerate(promoters):
            prom_file.write(f">promoter_{i}\n{record.seq[location.start:location.end]}\n")

    with open(non_promoter_output, "w") as non_prom_file:
        for i, location in enumerate(non_promoters):
            non_prom_file.write(f">non_promoter_{i}\n{record.seq[location.start:location.end]}\n")

    print(f"Saved {len(promoters)} promoters and {len(non_promoters)} non-promoters.")

find_promotors_with_parser("KU710884.1.gff", "KU710884.1.fasta", "result_folder")