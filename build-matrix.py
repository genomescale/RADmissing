import numpy
import os
import string
import random

# map specimen names to population names

specimen_map = {
	"cydno_ERS074404": "Heliconius_cydno_cordula",
	"cydno_ERS074405": "Heliconius_cydno_cordula",
	"cydno_ERS074406": "Heliconius_cydno_cordula",
	"cydno_ERS074407": "Heliconius_cydno_cordula",
	"cydno_ERS074408": "Heliconius_cydno_cordula",
	"cydno_ERS074409": "Heliconius_cydno_cordula",
	"elevatus_ERS070236": "Heliconius_elevatus",
	"elevatus_ERS070237": "Heliconius_elevatus",
	"elevatus_ERS070238": "Heliconius_elevatus",
	"elevatus_ERS070239": "Heliconius_elevatus",
	"elevatus_ERS070240": "Heliconius_elevatus",
	"ethilla_ERS070241": "Heliconius_ethilla_aerotome",
	"ethilla_ERS070242": "Heliconius_ethilla_aerotome",
	"ethilla_ERS070243": "Heliconius_ethilla_aerotome",
	"ethilla_ERS070244": "Heliconius_ethilla_aerotome",
	"heurippa_ERS074398": "Heliconius_heurippa",
	"heurippa_ERS074399": "Heliconius_heurippa",
	"heurippa_ERS074400": "Heliconius_heurippa",
	"heurippa_ERS074401": "Heliconius_heurippa",
	"heurippa_ERS074402": "Heliconius_heurippa",
	"heurippa_ERS074403": "Heliconius_heurippa",
	"melpomene_ERS070268": "Heliconius_melpomene_maletti",
	"melpomene_ERS070269": "Heliconius_melpomene_maletti",
	"melpomene_ERS070270": "Heliconius_melpomene_maletti",
	"melpomene_ERS070271": "Heliconius_melpomene_maletti",
	"melpomene_ERS070272": "Heliconius_melpomene_maletti",
	"melpomene_ERS070273": "Heliconius_melpomene_amaryllis",
	"melpomene_ERS070274": "Heliconius_melpomene_amaryllis",
	"melpomene_ERS070275": "Heliconius_melpomene_amaryllis",
	"melpomene_ERS070276": "Heliconius_melpomene_amaryllis",
	"melpomene_ERS070277": "Heliconius_melpomene_amaryllis",
	"pardalinus_ERS070257": "Heliconius_pardalinus_sergestus",
	"pardalinus_ERS070258": "Heliconius_pardalinus_sergestus",
	"pardalinus_ERS070259": "Heliconius_pardalinus_sergestus",
	"pardalinus_ERS070260": "Heliconius_pardalinus_sergestus",
	"pardalinus_ERS070261": "Heliconius_pardalinus_sergestus",
	"timareta_ERS070246": "Heliconius_timareta_peru",
	"timareta_ERS070247": "Heliconius_timareta_peru",
	"timareta_ERS070248": "Heliconius_timareta_peru",
	"timareta_ERS070249": "Heliconius_timareta_peru",
}

# sort specimens and populations into alphabetical order (separately).
# Later the alignment file is read as raw bytes, so we need to encode
# the specimen name strings as bytes first.

specimen_names = {specimen.encode() for specimen in specimen_map}
specimen_order = sorted(specimen_names)
population_order = sorted(set(specimen_map.values()))

n_specimens = len(specimen_order)
n_populations = len(population_order)

# build a map of specimen index to population index based on the previously
# constructed orders of specimens and populations

int_specimen_map = numpy.zeros(n_specimens, dtype = numpy.uint8)

for specimen_i, specimen_name in enumerate(specimen_order):
	population = specimen_map[specimen_name.decode()]
	population_i = population_order.index(population)
	int_specimen_map[specimen_i] = population_i

# these are the IUPAC ambiguity codes. Any other codes will be treated as
# missing data (e.g. Ns or dashes)

ambiguity_codes = {
	"A": ("A", "A"),
	"C": ("C", "C"),
	"G": ("G", "G"),
	"T": ("T", "T"),
	"R": ("A", "G"),
	"Y": ("C", "T"),
	"M": ("A", "C"),
	"K": ("G", "T"),
	"S": ("C", "G"),
	"W": ("A", "T"),
}

# build a map of ASCII codepoints for ambiguity codes to specific nucleotide
# letters. Uppercase letters begin at ASCII code point 65, so we need to
# offset the index of a letter in the alphabet (e.g. 0 for A) by this value
# (e.g. the ASCII codepoint for A is 0 + 65 = 65).

ascii_numeric_offset = 48
ascii_uppercase_offset = 65

int_ambiguity_codes = {}
for amb_code in ambiguity_codes:
	nt_code_a, nt_code_b = ambiguity_codes[amb_code]

	int_amb_code = string.ascii_uppercase.find(amb_code) + ascii_uppercase_offset
	int_nt_code_a = string.ascii_uppercase.find(nt_code_a) + ascii_uppercase_offset
	int_nt_code_b = string.ascii_uppercase.find(nt_code_b) + ascii_uppercase_offset

	int_ambiguity_codes[int_amb_code] = (int_nt_code_a, int_nt_code_b)

int_nt_codes = {}
for nt_int, nt in enumerate("ACGT"):
	ascii_codepoint = string.ascii_uppercase.find(nt) + ascii_uppercase_offset
	int_nt_codes[ascii_codepoint] = nt_int

# we can specify an arbitrary number of specimens sampled for each population.
# Because each specimen is diploid, the number of randomly sampled nucleotides
# will be double this.

n_specimens_per_population = 2
n_sequences_per_population = n_specimens_per_population * 2
n_sampled_specimens = n_populations * n_specimens_per_population

phy_path = "empirical_5/outfiles/empirical_5_m8.phy"
phy_file = open(phy_path, "rb")

is_matrix = False
alignments = []

alignment_matrix = None

locus_i = 0
specimen_i = 0

phy_header = phy_file.readline().strip().split()
concatenated_length = int(phy_header[1])

concatenated_matrix = numpy.zeros((n_specimens, concatenated_length), dtype = numpy.uint8)

l = phy_file.readline()
while l != b"":
	specimen_name, sequence = l.strip().split()
	if specimen_name in specimen_names:
		specimen_i = specimen_order.index(specimen_name)
		concatenated_matrix[specimen_i] = numpy.frombuffer(sequence, dtype = numpy.uint8)

	l = phy_file.readline()

partitions_path = "empirical_5/outfiles/empirical_5_m8.phy.partitions"
partitions_file = open(partitions_path)

# transposed matrix so the major dimension is columns and the minor dimension
# is specimens
concatenated_matrix_t = concatenated_matrix.transpose()

# read in and convert partition coordinates to C ranges (zero-indexing,
# half-open) which we call "partition boundaries"

partition_boundaries = []
for l in partitions_file.readlines():
	partition_coordinates = l.rstrip()[l.rfind("=") + 1:]
	partition_start, partition_end = partition_coordinates.split("-")
	c_partition_start = int(partition_start) - 1
	c_partition_end = int(partition_end)
	partition_boundaries.append((c_partition_start, c_partition_end))

n_partitions = len(partition_boundaries)

# a transposed matrix of unlinked SNPs with one row for each partition and one
# column for each sampled specimen. Where no biallelic columns have been
# sampled for a partition, that row will be all zeros. Where biallelic
# columns have been sampled, its index will be in included_partitions

snp_matrix_t = numpy.zeros((n_partitions, n_sampled_specimens), dtype = numpy.uint8)
included_partitions = []

# we can reuse these arrays rather than reallocating them

base_call_counts = numpy.zeros(n_populations, dtype = numpy.uint8)
population_offsets = numpy.zeros(n_populations, dtype = numpy.uint8)
sampled_column = numpy.zeros(n_sampled_specimens, dtype = numpy.uint8)
organized_column = numpy.zeros(n_specimens, dtype = numpy.uint8)

for partition_i, partition_coordinates in enumerate(partition_boundaries):
	c_partition_start, c_partition_end = partition_coordinates
	alignment_t = concatenated_matrix_t[c_partition_start:c_partition_end]

	# figure out if there are base calls for at least n specimens per population.
	# If there are, randomly sample n specimens pre population to create a
	# sampled column. Figure out if the resulting column is biallelic, and if
	# it is then add it to the list of biallelic columns

	biallelic_columns = []
	for column in alignment_t:
		base_call_counts.fill(0)
		population_offsets.fill(0)
		for specimen_i in range(n_specimens):
			population_i = int_specimen_map[specimen_i]
			if column[specimen_i] in int_ambiguity_codes:
				base_call_counts[population_i] += 1
				population_offsets[population_i + 1:n_populations] += 1

		if min(base_call_counts) >= n_specimens_per_population:
			for specimen_i in range(n_specimens):
				if column[specimen_i] in int_ambiguity_codes:
					population_i = int_specimen_map[specimen_i]
					population_i_offset = population_offsets[population_i]
					organized_column[population_i_offset] = column[specimen_i]
					population_offsets[population_i] += 1

			organized_population_range_start = 0
			sampled_population_range_start = 0
			sampled_nt_codes = set()
			for population_i in range(n_populations):
				organized_population_range_end = organized_population_range_start + base_call_counts[population_i]
				sampled_population_range_end = sampled_population_range_start + n_specimens_per_population

				population_ambiguity_codes = organized_column[organized_population_range_start:organized_population_range_end]
				sampled_ambiguity_codes = numpy.random.choice(population_ambiguity_codes, size = n_specimens_per_population, replace = False)
				sampled_column[sampled_population_range_start:sampled_population_range_end] = sampled_ambiguity_codes

				for int_amb_code in sampled_ambiguity_codes:
					sampled_nt_codes.update(int_ambiguity_codes[int_amb_code])

				organized_population_range_start = organized_population_range_end
				sampled_population_range_start = sampled_population_range_end

			if len(sampled_nt_codes) == 2:
				biallelic_columns.append(sampled_column)
				# allocate a new sampled column because we have to preserve the column
				# which is now in the list
				sampled_column = numpy.zeros(n_sampled_specimens, dtype = numpy.uint8)

	# if there are any biallelic columns, add one at random to the SNP matrix

	if len(biallelic_columns) > 0:
		random_biallelic_column = random.choice(biallelic_columns)
		snp_matrix_t[partition_i] = random_biallelic_column
		included_partitions.append(partition_i)

n_biallelic_markers = len(included_partitions)

# randomly resolve diploid codes into haploid bases, and record which bases
# exist for each biallelic column

binary_matrix_t = numpy.zeros((n_biallelic_markers, n_sampled_specimens * 2), dtype = numpy.uint8)

for biallelic_marker, partition_i in enumerate(included_partitions):
	snp_column = snp_matrix_t[partition_i]
	binary_column = binary_matrix_t[biallelic_marker]

	marker_nucleotides = set()

	sample_offset_start = 0
	for sample_i in range(n_sampled_specimens):
		sample_offset_end = sample_offset_start + 2
		sample_amb_code = snp_column[sample_i]
		sample_nt_codes = int_ambiguity_codes[sample_amb_code]
		marker_nucleotides.update(sample_nt_codes)
		binary_column[sample_offset_start:sample_offset_end] = sample_nt_codes
		sample_offset_start = sample_offset_end

	# recode bases as zeros and ones

	sorted_marker_nucleotides = sorted(marker_nucleotides)
	zero_base_i = numpy.random.randint(2)
	zero_base = sorted_marker_nucleotides[zero_base_i]

	for nucleotide_i, nucleotide in enumerate(binary_column):
		if nucleotide == zero_base:
			binary_column[nucleotide_i] = ascii_numeric_offset
		else:
			binary_column[nucleotide_i] = ascii_numeric_offset + 1

binary_matrix = binary_matrix_t.transpose()

binary_matrix_file = open("binary_matrix.fasta", "wb")

for sequence_i, sequence in enumerate(binary_matrix):
	population_i = sequence_i // (n_sequences_per_population)
	fake_taxon_i = sequence_i % (n_sequences_per_population)
	fake_taxon_letter_code = string.ascii_lowercase[fake_taxon_i]

	sequence_header = ">%s_%s\n" % (population_order[population_i], fake_taxon_letter_code)
	binary_matrix_file.write(sequence_header.encode())
	binary_matrix_file.write(sequence.tobytes())
	binary_matrix_file.write(b"\n")

binary_matrix_file.close()
