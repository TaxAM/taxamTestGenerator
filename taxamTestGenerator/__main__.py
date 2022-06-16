import numpy as np
import random
import csv
import os
import glob 
import sys
from utils.local_parsers import local_parsers

args = local_parsers().__dict__.copy()

random.seed(0)

def name_taxon(t, level=None):
  prefixs = ["RE", "FI", "CL", "OR", "FA", "GE", "ES"]
  if level is None:
    return ["NA" if t[l] == -1 else prefixs[l] + str(t[l]+1) for l in range(7)]
  return "NA" if t == -1 else prefixs[level] + str(t+1)

def generate_sample(
  sample_name,
  n_taxa_per_level,
  perc_partial_taxa,
  n_reads,
  n_contigs,
  perc_mapped_reads,
  n_taxa_reads,
  n_taxa_contigs,
  perc_class_reads,
  perc_class_contigs,
  perc_matched_class
):
  # /////////////////////////////////////////////////////////////////////////////
  # INNER FUNCTIONS
  # /////////////////////////////////////////////////////////////////////////////

  # Get the value of "name" in local and check if it is from 0 to 1.
  def valid_perc(name, locals = locals()):
    """Check if the value of "name" in local is from 0 to 1.

    Parameters
    ----------
    name : str
        Name of the variable to check.
    locals : dict, optional
        Dictionary of local variables, by default locals()

    Returns
    -------
    bool
        True if the value of "name" in local is from 0 to 1.
    """    
    if 0 <= eval(name, locals) <= 1:
      return True
    print(name + " should be a real number ranging from 0 to 1")

  def generate_taxa(number_of_taxa):
    """Generate a list of taxa.

    Parameters
    ----------
    number_of_taxa : int
        Number of taxa tuples in taxa list

    Returns
    -------
    list
        taxa tuples like [(5, 2, 7, 1, 2, 8, 3), (3, 1, 8, 6, 6, 5, 6), (2, 6, 8, 6, 0, 3, 0)] 
    """    
    taxa = []
    qtd = 0
    while True:
      # Tuple with 7 elements like (2, 7, 3, 8, 2, 3, 3)
      taxon = tuple([random.randint(0, n_taxa_per_level[i]-1) for i in range(7)])
      if taxon not in taxa:
        taxa.append(taxon)
        qtd += 1
        if qtd >= number_of_taxa:
          break

    changed = []
    # NA
    for _ in range(int(perc_partial_taxa * number_of_taxa)):
      target = None
      while True:
        target = random.randint(0, number_of_taxa-1)
        if target not in changed:
          break
      t = taxa.pop(target)[:random.randint(1, 5)]
      t = tuple([t[i] if i < len(t) else -1 for i in range(7)])
      taxa.insert(target, t)
      changed.append(target)
    return taxa

  def generate_inputs():
    if n_taxa_reads > 0:
      with open("reads_" + sample_name + ".tsv", 'w', newline='', encoding='utf-8') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        for id, taxon in read_classifications.items():
          row = ['READ' + str(id+1)]
          row.extend(name_taxon(taxon))
          csv_writer.writerow(row)

    with open("contigs_" + sample_name + ".tsv", 'w', newline='', encoding='utf-8') as f:
      csv_writer = csv.writer(f, delimiter='\t')
      for id, taxon in contig_classifications.items():
        row = ['CONTIG' + str(id+1)]
        row.extend(name_taxon(taxon))
        csv_writer.writerow(row)

    with open("mapping_" + sample_name + ".tsv", 'w', newline='', encoding='utf-8') as f:
      csv_writer = csv.writer(f, delimiter='\t')
      for r, c in mappings.items():
        row = ['READ' + str(r+1), 'CONTIG' + str(c+1)]
        csv_writer.writerow(row)

  def generate_outputs():
    for tie_strategy in ["read", "contig", "discard"]:
      for i in range(7):
        generate_output(i, tie_strategy)

    generate_inputs()

  def generate_output(level, tie_strategy):
    counter = {}
    for read_id in read_ids:
      taxon_read = read_classifications.get(read_id, [-1]*7)[:level+1]
      taxon_contig = contig_classifications.get(mappings.get(read_id, -1), [-1]*7)[:level+1]
      if taxon_contig[-1] < 0 and taxon_read[-1] < 0:
        continue
      if taxon_read[-1] < 0 or taxon_contig == taxon_read:
        counter[taxon_contig] = counter.get(taxon_contig, 0) + 1
      elif taxon_contig[-1] < 0:
        counter[taxon_read] = counter.get(taxon_read, 0) + 1
      else:
        if tie_strategy == "read":
          counter[taxon_read] = counter.get(taxon_read, 0) + 1
        if tie_strategy == "contig":
          counter[taxon_contig] = counter.get(taxon_contig, 0) + 1
    ordered_taxa = list(counter.keys())
    ordered_taxa.sort()
    with open("output_" + sample_name + "_level" + str(level) + "_tie-" + tie_strategy + ".tsv", 'w', newline='', encoding='utf-8') as f:
      csv_writer = csv.writer(f, delimiter='\t')
      for t in ordered_taxa:
        q = counter[t]
        csv_writer.writerow([t[-1], q])
    
  # /////////////////////////////////////////////////////////////////////////////
  # MAIN
  # /////////////////////////////////////////////////////////////////////////////
  # continue

  if len(n_taxa_per_level) != 7 or sum(list(map(lambda x : 0 if str(x).isdigit() else 1, n_taxa_per_level))) != 0:
    print("n_taxa_per_level should be a list with 7 integer elements")
    return -1
  n_taxa_per_level = list(map(lambda x : int(x), n_taxa_per_level))
  max_taxa = np.prod(n_taxa_per_level)
  if not valid_perc("perc_partial_taxa"):
    return -2
  if not valid_perc("perc_mapped_reads"):
    return -3
  if not valid_perc("perc_class_reads"):
    return -4
  if not valid_perc("perc_class_contigs"):
    return -5
  if not valid_perc("perc_matched_class"):
    return -6
  if n_taxa_reads > max_taxa:
    print("n_taxa_reads should not exceed the maximum number of taxa")
    return -7
  if n_taxa_contigs > max_taxa:
    print("n_taxa_reads should not exceed the maximum number of taxa")
    return -7
  if n_reads <= 0:
    print("n_reads should be a positive integer")
    return -8
  if n_contigs <= 0:
    print("n_contigs should be a positive integer")
    return -9
  read_ids = list(range(n_reads))
  contig_ids = list(range(n_contigs))
  mapped_read_ids = random.sample(read_ids, int(perc_mapped_reads * n_reads))
  mappings = {}
  
  for r in mapped_read_ids:
    mappings[r] = random.sample(contig_ids, 1)[0]
  all_taxa = generate_taxa(max(n_taxa_reads, n_taxa_contigs))

  taxa_contigs = random.sample(all_taxa, n_taxa_contigs)
  class_contig_ids = random.sample(read_ids, int(perc_class_contigs * n_reads))
  contig_classifications = {}
  for contig_id in random.sample(contig_ids, int(perc_class_contigs * n_contigs)):
    contig_classifications[contig_id] = random.sample(taxa_contigs, 1)[0]

  matching_buffer = []
  taxa_reads = random.sample(all_taxa, n_taxa_reads)
  read_classifications = {}
  for read_id in random.sample(read_ids, int(perc_class_reads * n_reads)):
    read_classifications[read_id] = random.sample(taxa_reads, 1)[0]
    if read_id in mappings:
      if mappings[read_id] in contig_classifications:
        matching_buffer.append((read_id, contig_classifications[mappings[read_id]]))

  for read_id, taxon in random.sample(matching_buffer, int(len(matching_buffer) * perc_matched_class)):
    read_classifications[read_id] = taxon
  generate_outputs()
  
def join_samples(pool_name, sample_names, remove_temporary=True):
  def load_sample(name):
    result = {}
    file_name = "output_" + name + "_level" + str(level) + "_tie-" + tie_strategy + ".tsv"
    with open(file_name, 'r') as f:
      csv_reader = csv.reader(f, delimiter='\t')
      for row in csv_reader:
        result[int(row[0])] = int(row[1])
    if remove_temporary:
        os.remove(file_name)
    return result
  os.mkdir(pool_name)
  for level in range(7):
    for tie_strategy in ["read", "contig", "discard"]:
      counter = {}
      tax_names = set()
      for sample_name in sample_names:
        counter[sample_name] = load_sample(sample_name)
        tax_names = tax_names.union(counter[sample_name].keys())
      tax_names = list(tax_names)
      tax_names.sort()
      with open("output_level" + str(level) + "_tie-" + tie_strategy + ".tsv", 'w', newline='', encoding='utf-8') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        # I have modified first field to be equal to TaxAM
        row = ['TaxAM']
        row.extend([s for s in sample_names])
        csv_writer.writerow(row)
        for tax_name in tax_names:
          row = [name_taxon(tax_name, level)]
          row.extend([counter[s].get(tax_name, 0) for s in sample_names])
          csv_writer.writerow(row)
  for f in glob.glob("*.tsv"):
    os.rename(f, os.path.join(pool_name, f))

def generate_samples(
  pool_name,
  sample_names,
  n_taxa_per_level,
  perc_partial_taxa,
  n_reads,
  n_contigs,
  perc_mapped_reads,
  n_taxa_reads,
  n_taxa_contigs,
  perc_class_reads,
  perc_class_contigs,
  perc_matched_class
):
  # Check if a folder for this sample already exists
  if os.path.exists(pool_name):
    print("Output folder already exists.")
    sys.exit(1)
  
  # Find all files with "*.tsv"
  if len(list(glob.glob("*.tsv"))):
    print("Current folder already has TSV file(s).")
    sys.exit(1)
  
  for sample_name in sample_names:
    generate_sample(
      sample_name,
      n_taxa_per_level,
      perc_partial_taxa,
      n_reads,
      n_contigs,
      perc_mapped_reads,
      n_taxa_reads,
      n_taxa_contigs,
      perc_class_reads,
      perc_class_contigs,
      perc_matched_class
    )

def adjust_classifications(pool_name, pattern):
    for p in glob.glob(os.path.join(pool_name, pattern)):
        buffer = []
        with open(p) as f:
            csv_reader = csv.reader(f, delimiter='\t')
            for row in csv_reader:
                if len(row)>0:
                    buffer.append([row[0], "; ".join(row[1:])])
        old = p
        p = p.rsplit('.', 1)[0] + '.txt'
        with open(p, 'w', newline='', encoding='utf-8') as f:
            csv_writer = csv.writer(f, delimiter='\t')
            for row in buffer:
                csv_writer.writerow(["C", row[0], 0, row[1]])
        os.remove(old)

# if len(sys.argv) != 13:
#     print("You must inform 12 parameters.")
#     sys.exit(1)

# Test
# python taxamTestGenerator pool_esc_a A,B 9,9,9,9,9,9,9 0 100 100 0.85 3000 1000 0.75 0.90 0.65
# python taxamTestGenerator -n pool_esc_a -s A,B -nt 9,9,9,9,9,9,9 -pt 0 -nr 100 -nc 100 -pm 0.85 -tr 3000 -tc 1000 -cr 0.75 -cc 0.90 -mc 0.65

# Inputs
pn = args['pool_name']
sn = args['sample_names'].split(",")
level = [int(x) for x in args['number_of_taxa_per_level'].split(",")]
generate_samples(
  pool_name = pn,
  sample_names = sn, 
  n_taxa_per_level = level,
  perc_partial_taxa = args['percentage_of_partial_taxa'],
  n_reads = args['number_of_reads'], 
  n_contigs = args['number_of_contigs'],
  perc_mapped_reads = args['percentage_of_mapped_reads'],
  n_taxa_reads = args['number_of_taxa_per_read'],
  n_taxa_contigs = args['number_of_taxa_per_contig'],
  perc_class_reads = args['percentage_of_classified_reads'], 
  perc_class_contigs = args['percentage_of_classified_contigs'],
  perc_matched_class = args['percentage_of_matched_class']
  )

join_samples(pn, sn)
adjust_classifications(pn, 'reads_*.tsv')
adjust_classifications(pn, 'contigs_*.tsv')