__author__ = 'hofmann'

import os
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.Validator.validator import Validator


class TaxonomicProfile(Validator):
	def __init__(self, taxonomy, logfile=None, verbose=True, debug=False):
		super(TaxonomicProfile, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
		assert isinstance(taxonomy, NcbiTaxonomy)
		self._taxonomy = taxonomy
		self._filename_taxonomic_profile = "taxonomic_profile_{sample_index}.txt"

	def write_taxonomic_profile_from_abundance_files(
		self, metadata_table_all, list_of_file_paths, directory_output, sample_id=""):
		metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		for index_abundance, file_path in enumerate(list_of_file_paths):
			community_abundance = metadata_table.parse_file(file_path, column_names=False)
			file_path_output = os.path.join(directory_output, self._filename_taxonomic_profile.format(index_abundance))
			with open(file_path_output, 'w') as stream_output:
				self.write_taxonomic_profile(
					community_abundance,
					stream_output,
					metadata_table_all,
					sample_id)

	def write_taxonomic_profile(self, community_abundance, stream_output, metadata_table, sample_id=""):
		assert isinstance(metadata_table, MetadataTable)
		genome_abundance = {}
		total_abundance = 0.0

		# for community in community_abundance:
		# 	all_communities += community

		for genome_id, file_name, distribution, total_length in community_abundance:
			if genome_id in genome_abundance:
				raise "genome id '{}' is not unique!".format(genome_id)
			genome_abundance[genome_id] = float(distribution)*float(total_length)
			total_abundance += genome_abundance[genome_id]

		for key, value in genome_abundance.items():
			genome_abundance[key] = value / total_abundance

			self._stream_taxonomic_profile(stream_output, genome_abundance, metadata_table, sample_id)

	def _stream_taxonomic_profile(self, stream_output, genome_id_to_percent, metadata_table, sample_id=""):
		strain_id_to_genome_id = {}
		genome_id_to_strain_id = {}

		genome_id_to_taxid = metadata_table.get_map(key_header="genome_ID", value_header="NCBI_ID")
		genome_id_to_otu = metadata_table.get_map(key_header="genome_ID", value_header="OTU")

		column_genome_id = metadata_table.get_column("genome_ID")
		column_strain_id = metadata_table.get_column("strain_id")
		if column_strain_id is None:
			column_strain_id = metadata_table.get_empty_column()
		else:
			genome_id_to_strain_id = metadata_table.get_map(key_header="genome_ID", value_header="strain_id")

		genome_id_to_lineage = self._get_genome_id_to_lineage(
			genome_id_to_percent.keys(), genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id)

		percent_by_rank_by_taxid = self._get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent)

		# add strain_id to metadata
		for row_index, genome_id in enumerate(column_genome_id):
			column_strain_id[row_index] = genome_id_to_strain_id[genome_id]
		assert len(column_strain_id) == len(set(column_strain_id))
		metadata_table.insert_column(column_strain_id, "strain_id")

		# stream taxonomic profile
		self._stream_tp_header(stream_output, sample_id)
		self._stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu)

	def _get_genome_id_to_lineage(
		self, list_of_genome_id, genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id):
		strains_by_taxid = {}
		genome_id_to_lineage = {}
		for genome_id in list_of_genome_id:
			tax_id = genome_id_to_taxid[genome_id]
			if tax_id == "":
				raise KeyError("genome_ID '{}' has no taxid!".format(genome_id))
			tax_id = self._taxonomy.get_updated_taxid(tax_id)
			genome_id_to_lineage[genome_id] = self._taxonomy.get_lineage_of_legal_ranks(
				tax_id, ranks=self._ranks, default_value=None)
			if genome_id_to_lineage[genome_id][-1] is not None:
				continue
			if tax_id not in strains_by_taxid:
				strains_by_taxid[tax_id] = 0
			strains_by_taxid[tax_id] += 1
			if genome_id in genome_id_to_strain_id and genome_id_to_strain_id[genome_id]:
				strain_id = genome_id_to_strain_id[genome_id]
			else:
				strain_id = "{}.{}".format(tax_id, strains_by_taxid[tax_id])
				genome_id_to_strain_id[genome_id] = strain_id
			genome_id_to_lineage[genome_id][-1] = strain_id
			strain_id_to_genome_id[strain_id] = genome_id
		return genome_id_to_lineage

	def _get_percent_by_rank_by_taxid(self, genome_id_to_lineage, genome_id_to_percent):
		percent_by_rank_by_taxid = {}
		for rank in self._ranks:
			percent_by_rank_by_taxid[rank] = dict()

		for rank_index, rank in enumerate(self._ranks):
			# rank = ranks[rank_index]
			for genome_id in genome_id_to_lineage:
				tax_id = genome_id_to_lineage[genome_id][rank_index]
				if tax_id is None:
					continue
				percent = genome_id_to_percent[genome_id]
				if tax_id not in percent_by_rank_by_taxid[rank]:
					percent_by_rank_by_taxid[rank][tax_id] = 0
				percent_by_rank_by_taxid[rank][tax_id] += percent
		return percent_by_rank_by_taxid

	def _stream_tp_rows(self, stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu):
		row_format = "{taxid}\t{rank}\t{taxpath}\t{taxpath_sn}\t{abp:.4f}\t{gid}\t{otu}\n"
		for rank_index, rank in enumerate(self._ranks):
			for tax_id in percent_by_rank_by_taxid[rank]:
				if '.' in tax_id:
					genome_id = strain_id_to_genome_id[tax_id]
					otu = genome_id_to_otu[genome_id]
					lineage = self._taxonomy.get_lineage_of_legal_ranks(tax_id.split('.')[0], ranks=self._ranks, default_value="")
					lineage[-1] = tax_id
				else:
					genome_id = ""
					otu = ""
					lineage = self._taxonomy.get_lineage_of_legal_ranks(tax_id, ranks=self._ranks, default_value="")

				lineage = lineage[:rank_index+1]
				lineage_sn = [self._taxonomy.get_scientific_name(tid) if tid != "" and '.' not in tid else "" for tid in lineage]
				if '.' in tax_id:
					lineage_sn[-1] = self._taxonomy.get_scientific_name(tax_id.split('.')[0]) + " strain"  # ""

				stream_output.write(row_format.format(
					taxid=tax_id,
					rank=rank,
					taxpath="|".join(lineage),
					taxpath_sn="|".join(lineage_sn),
					abp=percent_by_rank_by_taxid[rank][tax_id]*100,
					gid=genome_id,
					otu=otu
				))

	def _stream_tp_header(self, output_stream, identifier):
		output_stream.write("@SampleID:{}\n".format(identifier))
		output_stream.write("@Version:0.9.1\n")
		output_stream.write("@Ranks:{ranks}\n\n".format(ranks="|".join(self._ranks)))
		output_stream.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\t_CAMI_genomeID\t_CAMI_OTU\n")

	@staticmethod
	def print_ratios(communities):
		com_length = [0] * len(communities)
		total_length = 0
		for index, community in enumerate(communities):
			for line in community:
				com_length[index] += float(line[3]) * float(line[4])
			total_length += com_length[index]
		print("community length:", com_length)
		for x in com_length:
			print("community ratio:", x/total_length)
