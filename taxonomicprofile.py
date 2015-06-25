__author__ = 'hofmann'

import os
from communitydesign import Community
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy


class TaxonomicProfile(object):
	def __init__(self, taxonomy):
		self._ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
		assert isinstance(taxonomy, NcbiTaxonomy)
		self._taxonomy = taxonomy

	def main(self, list_of_communities, directory_distributions, directory_output, number_of_samples, metadata_table_all):
		"""

		@param list_of_communities:
		@type list_of_communities: list[Community]
		@param directory_distributions:
		@type :
		@param metadata_table_all:
		@type :
		@return:
		@rtype:
		"""

		assert isinstance(list_of_communities, list)
		for community in list_of_communities:
			assert isinstance(community, Community)
		assert isinstance(metadata_table_all, MetadataTable)
		file_path_prefix = os.path.join(directory_distributions, options.file_prefix)
		taxonomic_file_path_prefix = os.path.join(directory_output, options.file_prefix)

		# read communities and adapt to ratio
		for j in range(number_of_samples):  # TODO sample number
			communities = []
			communities_length = []
			total_length = 0
			for i in range(len(list_of_communities)):
				dist_file_name = file_path_prefix+'_'+str(i)+'_'+str(j)+'.txt'
				communities.append(self.read_community(dist_file_name))
				communities_length.append(0)
				os.remove(dist_file_name)

				genomes = set()
				for sequence_info in communities[i]:
					if sequence_info[0] in genomes:
						continue
					genomes.add(sequence_info[0])
					communities_length[i] += float(sequence_info[3])  # * float(sequence_info[4])
				total_length += communities_length[i]

			# out.append(read_communities[0][0])
			for i in range(len(list_of_communities)):
				# factor = float((float(options.communities[i]['ratio'])/float(options.communities[0]['ratio']))*
				# communities[0][1]/communities[i][1])#new_sum/read_communities[i][1]
				# factor = ratio / (comunity_length / total_length)
				factor = float(list_of_communities[i]['ratio']) / (communities_length[i] / float(total_length))
				self.update_community(communities[i], factor)

			# print_ratios(communities)
			self.save_community(communities, file_path_prefix+'_'+str(j)+'.txt')
			sample_id = ""  # "{}{}".format(options.dataset_id, j)

			if self._taxonomy is not None:
				self.save_community_profile(
					communities, taxonomic_file_path_prefix+'_taxonomic_profile_'+str(j)+'.tsv',
					metadata_table_all, sample_id)
		return True

	@staticmethod
	def read_community(file_path):
		lines = []
		try:
			# print(file_path)
			with open(file_path, 'r') as community_file:
				for line in community_file:
					lines.append(line.strip().split('\t'))
		except:
			raise Exception
			# lines[-1].append(float(lines[-1][3]) * float(lines[-1][4]))
		return lines

	@staticmethod
	def update_community(community, factor):
		for index in range(len(community)):
			# new_abundance = old_abundance * factor
			community[index][3] = str(float(community[index][3])*factor)

	@staticmethod
	def save_community(communities, path):
		try:
			with open(path, 'w') as new_file:
				for community in communities:
					for line in community:
						new_file.write('\t'.join(line) + '\n')
		except Exception as e:
			raise e

	def save_community_profile(self, communities, path, metatable, sample_id=""):
		assert isinstance(metatable, MetadataTable)
		all_communities = []
		genome_percentage = {}
		percentage_sum = 0.0

		for community in communities:
			all_communities += community

		for entry in all_communities:
			genome_id, sequence_id, file_name, abundance, sequence_length, dna_type = entry

			tmp = float(abundance)*float(sequence_length)
			if genome_id in genome_percentage:
				genome_percentage[genome_id] += tmp
			else:
				genome_percentage[genome_id] = tmp

			percentage_sum += tmp

		for key, value in genome_percentage.items():
			genome_percentage[key] = value / percentage_sum

		# try:
		with open(path, 'w') as new_file:
			self.write_taxonomic_profile(new_file, genome_percentage, metatable, sample_id)
		# except Exception as e:
		#    raise e

	def write_taxonomic_profile(self, stream_output, genome_id_to_percent, metadata_table, sample_id):
		strain_id_to_genome_id = {}
		genome_id_to_strain_id = {}

		genome_id_to_taxid = metadata_table.get_map(key_header="genome_ID", value_header="NCBI_ID")
		genome_id_to_otu = metadata_table.get_map(key_header="genome_ID", value_header="OTU")

		# add strain_id to metadata
		column_genome_id = metadata_table.get_column("genome_ID")
		column_strain_id = metadata_table.get_empty_column()
		for row_index in range(len(column_genome_id)):
			genome_id = column_genome_id[row_index]
			column_strain_id[row_index] = genome_id_to_strain_id[genome_id]
		metadata_table.insert_column(column_strain_id, "strain_id")

		genome_id_to_lineage = self._get_genome_id_to_lineage(
			genome_id_to_percent.keys(), genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id)

		percent_by_rank_by_taxid = self._get_percent_by_rank_by_taxid(genome_id_to_lineage, genome_id_to_percent)

		self._stream_tp_header(stream_output, sample_id)
		self._stream_tp_rows(stream_output, percent_by_rank_by_taxid, strain_id_to_genome_id, genome_id_to_otu)

	def _get_genome_id_to_lineage(
		self, list_of_genome_id, genome_id_to_taxid, strain_id_to_genome_id, genome_id_to_strain_id):
		strains_by_taxid = {}
		genome_id_to_lineage = {}
		for genome_id in list_of_genome_id:
			taxid = genome_id_to_taxid[genome_id]
			if taxid == "":
				raise KeyError("genome_ID '{}' has no taxid!".format(genome_id))
			genome_id_to_lineage[genome_id] = self._taxonomy.get_lineage_of_legal_ranks(
				taxid, ranks=self._ranks, default_value=None)
			if genome_id_to_lineage[genome_id][-1] is not None:
				continue
			if taxid not in strains_by_taxid:
				strains_by_taxid[taxid] = 0
			strains_by_taxid[taxid] += 1
			strain_id = "{}.{}".format(taxid, strains_by_taxid[taxid])
			genome_id_to_lineage[genome_id][-1] = strain_id
			strain_id_to_genome_id[strain_id] = genome_id
			genome_id_to_strain_id[genome_id] = strain_id
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
					lineage[-1] = tax_id  # ""
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
