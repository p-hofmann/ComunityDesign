__author__ = 'hofmann'


import os
import random
import tempfile
from Bio import SeqIO
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.Validator.sequencevalidator import SequenceValidator
from scripts.StrainSimulationWrapper.strainsimulationwrapper import StrainSimulationWrapper
from scripts.StrainSelector.strainselector import StrainSelector
from scripts.PopulationDistribution.populationdistribution import PopulationDistribution


# ##################################
#
#          Community
#
# ##################################


class Community(object):

	def __init__(
		self, identifier, genomes_total, genomes_real, limit_per_otu, file_path_metadata_table,
		file_path_genome_locations, file_path_gff_locations, ratio, mode,
		log_mu, log_sigma, gauss_mu=None, gauss_sigma=None,
		verbose=False):
		"""

		@param identifier:
		@type identifier: str | unicode
		@param genomes_total:
		@type genomes_total: int
		@param genomes_real:
		@type genomes_real: int
		@param limit_per_otu:
		@type limit_per_otu: int
		@param file_path_metadata_table:
		@type file_path_metadata_table: str | unicode
		@param file_path_genome_locations:
		@type file_path_genome_locations: str | unicode
		@param file_path_gff_locations:
		@type file_path_gff_locations: str | unicode
		@param ratio:
		@type ratio: int | long | float
		@param mode:
		@type mode: str | unicode
		@param log_mu:
		@type log_mu: int | long | float
		@param log_sigma:
		@type log_sigma: int | long | float
		@param gauss_mu:
		@type gauss_mu: int | long | float
		@param gauss_sigma:
		@type gauss_sigma: int | long | float
		@param verbose:
		@type verbose: bool
		"""
		assert genomes_real is None or genomes_real <= genomes_total
		self.genomes_real = genomes_real
		self.genomes_total = genomes_total
		self.limit_per_otu = limit_per_otu
		self.file_path_metadata_table = file_path_metadata_table
		self.file_path_genome_locations = file_path_genome_locations
		self.file_path_gff_locations = file_path_gff_locations
		self.ratio = ratio
		self.log_mu = log_mu
		self.log_sigma = log_sigma
		self.gauss_mu = gauss_mu
		self.gauss_sigma = gauss_sigma
		self.mode = mode
		self.simulate_strains = False
		if genomes_real and genomes_real <= genomes_total:
			self.simulate_strains = True
		self.verbose = verbose
		self.id = identifier


# ##################################
#
#          PrepareStrains
#
# ##################################


class PrepareStrains(SequenceValidator):
	_label = "PrepareStrains"

	def _get_genome_id_to_path_map(self, file_path_of_file_mapping_genome_id_to_paths, list_of_drawn_genome_id):
		"""

		@param file_path_of_file_mapping_genome_id_to_paths:
		@type file_path_of_file_mapping_genome_id_to_paths: str | unicode
		@param list_of_drawn_genome_id:
		@type list_of_drawn_genome_id: list[str|unicode]

		@return:
		@rtype: dict[str|unicode, str|unicode]
		"""
		mdt = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		mdt.read(file_path_of_file_mapping_genome_id_to_paths)
		genome_id_to_path_map = mdt.get_map(0, 1, unique_key=True)
		assert set(genome_id_to_path_map.keys()).issuperset(list_of_drawn_genome_id)
		return {genome_id: genome_id_to_path_map[genome_id] for genome_id in list_of_drawn_genome_id}

	def _move_genome_file(self, input_path, output_path, min_length, list_of_sequence_names=None, file_format="fasta"):
		"""

		@param input_path:
		@type input_path: str | unicode
		@param output_path:
		@type output_path: str | unicode
		@param min_length:
		@type min_length: int | long
		@param list_of_sequence_names:
		@type list_of_sequence_names: set[str|unicode]
		@param file_format:
		@type file_format: str | unicode

		@return: @raise Exception:
		@rtype: int | long
		"""
		assert self.validate_file(input_path)
		assert not self.validate_file(output_path, silent=True), "Overwriting files prohibited: '{}'".format(output_path)
		if list_of_sequence_names is None:
			list_of_sequence_names = []
		with open(input_path, 'r') as stream_input, open(output_path, 'w') as stream_output:
			total_base_pairs = self._cleanup_and_filter_sequences(
				stream_input, stream_output, min_length, list_of_sequence_names, file_format)
		if total_base_pairs == 0:
			msg = "No valid sequences in '{}'".format(stream_input.name)
			self._logger.error(msg)
			raise Exception(msg)
		return total_base_pairs

	def _cleanup_and_filter_sequences(
		self, stream_input, stream_output, min_length, list_of_sequence_names, file_format="fasta"):
		"""

		@param stream_input:
		@type stream_input: file | FileIO | StringIO
		@param stream_output:
		@type stream_output: file | FileIO | StringIO
		@param min_length:
		@type min_length: int | long
		@param list_of_sequence_names:
		@type list_of_sequence_names: set[str|unicode]
		@param file_format:
		@type file_format: str | unicode
		@return:
		@rtype: int | long
		"""
		total_base_pairs = 0
		for seq_record in SeqIO.parse(stream_input, file_format):
			# remove description, else art illumina messes up sam format
			seq_record.description = ''
			if len(seq_record.seq) < min_length:
				self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
					os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
				continue
			if seq_record.id in list_of_sequence_names:
				seq_record.id = self._get_new_name(seq_record.id, list_of_sequence_names)
			list_of_sequence_names.add(seq_record.id)
			# file_handler.write(">{}\n".format(sequence_id))
			# file_handler.writelines("{}\n".format(seq_record.seq))
			stream_output.write(seq_record.format(file_format))
			total_base_pairs += len(seq_record.seq)
		return total_base_pairs

	def move_genome_files(
		self, genome_id_to_path_map, directory_genomes, min_length, list_of_sequence_names=None):
		"""

		@param genome_id_to_path_map:
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]
		@param directory_genomes:
		@type directory_genomes: str | unicode
		@param min_length:
		@type min_length: int | long
		@param list_of_sequence_names:
		@type list_of_sequence_names: list[str|unicode]

		@return:
		@rtype: dict[str|unicode, int|long]
		"""
		assert isinstance(genome_id_to_path_map, dict)
		if list_of_sequence_names is None:
			list_of_sequence_names = set()
		genome_id_to_total_length = {}
		for genome_id, genome_file_path in genome_id_to_path_map.iteritems():
			file_name = os.path.basename(genome_file_path)
			new_genome_file_path = os.path.join(directory_genomes, file_name)
			total_length = self._move_genome_file(genome_file_path, new_genome_file_path, min_length, list_of_sequence_names)
			genome_id_to_total_length[genome_id] = total_length
		return genome_id_to_total_length

	@staticmethod
	def _get_new_name(name, list_of_sequence_names):
		"""

		@param name:
		@type name: str | unicode
		@param list_of_sequence_names:
		@type list_of_sequence_names: set[str|unicode]
		@return:
		@rtype : str | unicode
		"""
		index = 0
		new_name = name
		while new_name in list_of_sequence_names:
			new_name = name + "_{}".format(index)
			index += 1
		return new_name

	def validate_format(self, list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
		"""

		@param list_of_file_paths:
		@type list_of_file_paths: list[str|unicode]
		@param file_format:
		@type file_format: str | unicode
		@param sequence_type:
		@type sequence_type: str | unicode
		@param ambiguous:
		@type ambiguous: bool

		@return: True if all valid
		@rtype: bool
		"""
		result = True
		for file_path in list_of_file_paths:
			if not self.validate_sequence_file(file_path, file_format, sequence_type, ambiguous):
				result = False
		return result


# ##################################
#
#          CommunityDesign
#
# ##################################


class CommunityDesign(PrepareStrains):
	_label = "CommunityDesign"

	_filename_distribution_comunity = "distribution_{comunity_index}_{sample_index}.txt"
	_filename_distribution_comunity_joint = "distribution_{sample_index}.txt"

	# TODO: plasmids within genome files
	# used_genomes_with_plasmids[genome_id] = random.randint(7, 10)
	# distribution = str(int(distribution) * factor)
	def __init__(
		self, column_name_genome_id="genome_ID", column_name_otu="OTU",
		column_name_novelty_category="novelty_category", column_name_ncbi="NCBI_ID", column_name_source="source",
		max_processors=1, tmp_dir=None, logfile=None, verbose=True, debug=False, seed=None):
		"""

		@param column_name_genome_id:
		@type column_name_genome_id: str | unicode
		@param column_name_otu:
		@type column_name_otu: str | unicode
		@param column_name_novelty_category:
		@type column_name_novelty_category: str | unicode
		@param column_name_ncbi:
		@type column_name_ncbi: str | unicode
		@param column_name_source:
		@type column_name_source: str | unicode
		@param max_processors: maximum number of processors available to be used
		@type max_processors: long | int
		@param tmp_dir: working directory or place temporary files can be stored
		@type tmp_dir: str | unicode
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		super(CommunityDesign, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		if seed is not None:
			random.seed(seed)
		self._seed = seed
		# self._filename_distribution = filename_prefix_distribution + "{index}.txt"
		self._column_name_genome_id = column_name_genome_id
		self._column_name_otu = column_name_otu
		self._column_name_novelty_category = column_name_novelty_category
		self._column_name_source = column_name_source
		self._column_name_ncbi = column_name_ncbi

		assert isinstance(max_processors, (long, int))
		assert max_processors > 0
		self._max_processors = max_processors

		if tmp_dir is None:
			tmp_dir = tempfile.gettempdir()
		self._tmp_dir = tmp_dir
		assert self.validate_dir(self._tmp_dir)

	def _write_distribution_files(
		self, directory_output, genome_id_to_total_length, genome_id_to_distributions, genome_id_to_file_name,
		index_community):
		"""

		@param directory_output:
		@type directory_output: str | unicode
		@param genome_id_to_total_length:
		@type genome_id_to_total_length: dict[str|unicode, int|long]
		@param genome_id_to_distributions:
		@type genome_id_to_distributions: dict[str|unicode, list[float]]
		@param genome_id_to_file_name:
		@type genome_id_to_file_name: dict[str|unicode, str|unicode]
		@param index_community:
		@type index_community: int | str | unicode
		"""
		list_of_genome_id = genome_id_to_distributions.keys()
		distributions = len(genome_id_to_distributions[list_of_genome_id[0]])
		for index_sample in range(0, distributions):
			out_file_path = os.path.join(
				directory_output, self._filename_distribution_comunity.format(
					comunity_index=index_community, sample_index=index_sample))
			with open(out_file_path, 'w') as out_file_handler:
				for genome_id in genome_id_to_distributions:
					distribution = genome_id_to_distributions[genome_id][index_sample]
					file_name = genome_id_to_file_name[genome_id]
					total_length = genome_id_to_total_length[genome_id]
					out_file_handler.write("{id}\t{file_name}\t{dist}\t{len}\n".format(
						id=genome_id,
						file_name=file_name,
						dist=distribution,
						len=total_length))

	def design(
		self, community, number_of_samples, metadata_table_comby,
		directory_distributions, directory_genomes, directory_metadata, directory_template,
		list_of_sequence_names=None, min_length=1):
		"""

		@param community:
		@type community: Community
		@param number_of_samples:
		@type number_of_samples: int
		@param metadata_table_comby:
		@type metadata_table_comby: MetadataTable
		@param directory_distributions:
		@type directory_distributions: str | unicode
		@param directory_genomes:
		@type directory_genomes: str | unicode
		@param directory_metadata:
		@type directory_metadata: str | unicode
		@param directory_template:
		@type directory_template: str | unicode
		@param list_of_sequence_names:
		@type list_of_sequence_names: list[str|unicode]
		@param min_length:
		@type min_length: int | long

		@return:
		@rtype: list[str|unicode]
		"""
		assert isinstance(community, Community)
		assert isinstance(metadata_table_comby, MetadataTable)

		if list_of_sequence_names is None:
			list_of_sequence_names = set()

		number_of_strains = community.genomes_total

		# pick how much a strain will be simulated
		genome_amounts = []
		strain_simulation = None
		if community.simulate_strains:
			strain_simulation = StrainSimulationWrapper(
				executable_sim=None,
				directory_template=directory_template,
				column_name_gid=self._column_name_genome_id,
				column_name_ncbi=self._column_name_ncbi,
				column_name_source=self._column_name_source,
				separator='\t',
				filename_prefix="simulated_",
				keep_original=True,
				max_processors=self._max_processors,
				tmp_dir=self._tmp_dir,
				logfile=self._logfile, verbose=self._verbose, debug=self._debug, seed=self._seed)

			probability = None  # 1-options.communities[community_id]["evolve"]
			genome_amounts = strain_simulation.get_genome_amounts(
				probability=probability,
				max_genome_amount=community.genomes_total,
				num_real_genomes=community.genomes_real,
				silent=not community.verbose
			)
			number_of_strains = len(genome_amounts)

		# draw strains
		self._logger.info("Drawing strains.")
		metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		metadata_table.read(community.file_path_metadata_table, column_names=True)
		strain_selector = StrainSelector(
			column_name_genome_id=self._column_name_genome_id,
			column_name_otu=self._column_name_otu,
			column_name_novelty_category=self._column_name_novelty_category,
			logfile=self._logfile, verbose=self._verbose, debug=self._debug
			)
		list_of_drawn_genome_id = strain_selector.get_drawn_genome_id(
			metadata_table=metadata_table,
			number_of_strains=number_of_strains,
			number_of_strains_per_otu=community.limit_per_otu
			)

		# write unused data to separate file
		old_base_name = os.path.basename(community.file_path_metadata_table)
		file_prefix, extention = os.path.splitext(old_base_name)
		new_file_name = "unused_c{index}_{prefix}{ext}".format(
			prefix=file_prefix,
			index=community.id,
			ext=extention)
		metadata_new_file_path = os.path.join(directory_metadata, new_file_name)
		metadata_table.write(
			metadata_new_file_path,
			exclude=True,
			value_list=list_of_drawn_genome_id,
			key_column_name=self._column_name_genome_id,
			column_names=True)

		# get path for every genome
		genome_id_to_file_path_gff = self._get_genome_id_to_path_map(
			community.file_path_gff_locations, list_of_drawn_genome_id)
		genome_id_to_path_map = self._get_genome_id_to_path_map(
			community.file_path_genome_locations, list_of_drawn_genome_id)

		# concatenate
		metadata_table.reduce_rows_to_subset(list_of_drawn_genome_id, self._column_name_genome_id)
		metadata_table_comby.concatenate(metadata_table, strict=False)

		# validate correct format of files
		self._logger.info("Validating raw sequence files!")
		assert self.validate_format(
			list_of_file_paths=genome_id_to_path_map.values(),
			file_format="fasta",
			sequence_type="dna",
			ambiguous=True
			), "Validation of file format failed!"

		# simulate diversity around strains
		if community.simulate_strains:
			genome_id_to_amounts = strain_simulation.get_genome_id_to_amounts(list_of_drawn_genome_id, genome_amounts)
			strain_simulation.simulate_strains(
				meta_table=metadata_table_comby,
				genome_id_to_amounts=genome_id_to_amounts,
				genome_id_to_file_path_genome=genome_id_to_path_map,
				genome_id_to_file_path_gff=genome_id_to_file_path_gff)
			# adopt new list that includes simulated strains
			self._logger.info("Validating simulated sequence files!")
			for genome_id, file_path in genome_id_to_path_map.iteritems():
				if genome_id in list_of_drawn_genome_id:
					continue
				assert self.validate_sequence_file(
					file_path,
					file_format="fasta",
					sequence_type="dna",
					ambiguous=True)
			list_of_drawn_genome_id = genome_id_to_path_map.keys()

		# get community distributions
		population_distribution = PopulationDistribution(
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)
		list_of_distributions = population_distribution.get_lists_of_distributions(
			size_of_population=len(list_of_drawn_genome_id),
			number_of_samples=number_of_samples,
			modus=community.mode,
			log_mu=community.log_mu, log_sigma=community.log_sigma,
			gauss_mu=community.gauss_mu, gauss_sigma=community.gauss_sigma
		)

		# move and clean up files (removes sequence description)
		genome_id_to_total_length = self.move_genome_files(
			genome_id_to_path_map,
			directory_genomes=directory_genomes,
			min_length=min_length,
			list_of_sequence_names=list_of_sequence_names)

		# write distribution file
		genome_id_to_distributions = self._get_genome_id_to_distributions(list_of_drawn_genome_id, list_of_distributions)
		genome_id_to_file_name = self._get_genome_id_to_file_name(genome_id_to_path_map)
		self._write_distribution_files(
			directory_output=directory_distributions,
			genome_id_to_total_length=genome_id_to_total_length,
			genome_id_to_distributions=genome_id_to_distributions,
			genome_id_to_file_name=genome_id_to_file_name,
			index_community=community.id)
		return list_of_drawn_genome_id

	@staticmethod
	def _get_genome_id_to_distributions(list_of_drawn_genome_id, list_of_distributions):
		"""

		@param list_of_drawn_genome_id:
		@type list_of_drawn_genome_id: list[str|unicode]
		@param list_of_distributions:
		@type list_of_distributions: list[list[float]]

		@return:
		@rtype : dict[str|unicode, list[float]]
		"""
		genome_id_to_distributions = {}
		assert len(list_of_drawn_genome_id) == len(list_of_distributions)
		for index in range(len(list_of_drawn_genome_id)):
			genome_id_to_distributions[list_of_drawn_genome_id[index]] = list_of_distributions[index]
		return genome_id_to_distributions

	@staticmethod
	def _get_genome_id_to_file_name(genome_id_to_path_map):
		"""

		@param genome_id_to_path_map:
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]

		@return:
		@rtype : dict[str|unicode, str|unicode]
		"""
		genome_id_to_file_name = {}
		for genome_id, file_path in genome_id_to_path_map.iteritems():
			genome_id_to_file_name[genome_id] = os.path.basename(file_path)
		return genome_id_to_file_name

	#####################
	#
	# merge communities
	#
	#####################

	def merge_communities(self, list_of_communities, directory_distributions, number_of_samples, metadata_table_all):
		"""

		@param list_of_communities:
		@type list_of_communities: list[Community]
		@param directory_distributions:
		@type directory_distributions: str | unicode
		@param number_of_samples:
		@type number_of_samples: int | long
		@param metadata_table_all:
		@type metadata_table_all: MetadataTable

		@return:
		@rtype: list[str|unicode]
		"""
		assert isinstance(list_of_communities, list)
		for community in list_of_communities:
			assert isinstance(community, Community)
		assert isinstance(metadata_table_all, MetadataTable)

		# read communities and adapt to ratio
		list_of_output_paths = []
		for index_sample in range(number_of_samples):
			communities = []
			list_of_community_total_distribution = [0] * len(list_of_communities)
			sample_total_distribution = 0

			# create list of community files for one sample
			list_of_file_paths = []
			for index_community, c in enumerate(list_of_communities):
				dist_file_name = os.path.join(
					directory_distributions,
					self._filename_distribution_comunity.format(
						comunity_index=index_community, sample_index=index_sample))
				list_of_file_paths.append(dist_file_name)

			metadata_table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
			for index_community, file_path in enumerate(list_of_file_paths):
				communities.append(metadata_table.parse_file(file_path, column_names=False))
				# communities_length.append(0)

				genomes = set()
				for genome_id, file_name, distribution, genome_length in communities[index_community]:
					if genome_id in genomes:
						# print genome_id, filename, distribution, genome_length
						raise ValueError("Genome id '{}' not unique".format(genome_id))
					genomes.add(genome_id)
					list_of_community_total_distribution[index_community] += float(distribution)  # * float(sequence_info[4])
				sample_total_distribution += list_of_community_total_distribution[index_community]
				communities[index_community].close()

			# out.append(read_communities[0][0])
			list_of_community_factor = [0.0] * len(list_of_communities)
			for index_community, c in enumerate(list_of_file_paths):
				ratio = float(list_of_communities[index_community].ratio)
				community_total_distribution = float(list_of_community_total_distribution[index_community])
				current_proportion_in_sample = community_total_distribution / float(sample_total_distribution)
				list_of_community_factor[index_community] = ratio / current_proportion_in_sample
				# self.update_community(communities[index_community], factor)

			# join communities
			communities = []
			for index_community, file_path in enumerate(list_of_file_paths):
				communities.append(metadata_table.parse_file(file_path, column_names=False))

			# print_ratios(communities)
			file_path_output = os.path.join(
				directory_distributions, self._filename_distribution_comunity_joint.format(sample_index=index_sample))
			with open(file_path_output, 'w') as stream_output:
				self._write_joined_community(
					communities,
					list_of_community_factor,
					stream_output)
			list_of_output_paths.append(file_path_output)

			# delete now obsolete files
			for file_path in list_of_file_paths:
				os.remove(file_path)
		return list_of_output_paths

	@staticmethod
	def _write_joined_community(communities, list_of_community_factor, stream_output):
		"""

		@param communities:
		@type communities: generator[ dict[int|long|str|unicode, str|unicode] ]
		@param list_of_community_factor:
		@type list_of_community_factor: list[float]
		@param stream_output:
		@type stream_output: file | FileIO | StringIO | basestring
		"""
		line_format = "{gid}\t{filename}\t{distr}\t{length}\n"
		for community_index, community in enumerate(communities):
			factor = list_of_community_factor[community_index]
			for genome_id, filename, distribution, genome_length in community:
				stream_output.write(line_format.format(
					gid=genome_id,
					filename=filename,
					distr=float(distribution) * factor,
					length=genome_length,
				))
			community.close()
