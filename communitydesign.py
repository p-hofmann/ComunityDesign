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
		Accumulation of all community related information

		@param identifier: Community identifier
		@type identifier: str | unicode
		@param genomes_total: Total amount of genomes to be drawn from this community
		@type genomes_total: int
		@param genomes_real: Amount of real genomes to be drawn, rest will drawn from simulated ones
		@type genomes_real: int
		@param limit_per_otu: A Maximum for drawn genomes belonging to the same otu, unless more are required to be drawn
		@type limit_per_otu: int
		@param file_path_metadata_table: Table of Metadata for each genome of the community
		@type file_path_metadata_table: str | unicode
		@param file_path_genome_locations: Format: 'id \t file path to fasta file'
		@type file_path_genome_locations: str | unicode
		@param file_path_gff_locations: Format: 'id \t file path to gff file'
		@type file_path_gff_locations: str | unicode
		@param ratio: If one comm. has ratio=1 and another has ration=2, the other community will be twice the size
		@type ratio: int | long | float
		@param mode: Valid: 'replicates', 'timeseries_normal', 'timeseries_lognormal', 'differential'
		@type mode: str | unicode
		@param log_mu: Mean of drawn log distribution
		@type log_mu: int | long | float
		@param log_sigma: Standard deviation of log distribution
		@type log_sigma: int | long | float
		@param gauss_mu: Mean of drawn gauss distribution
		@type gauss_mu: int | long | float
		@param gauss_sigma: Standard deviation of gauss distribution
		@type gauss_sigma: int | long | float
		@param verbose: More output and user interaction is enabled.
		@type verbose: bool
		"""
		assert genomes_real is None or genomes_real <= genomes_total
		assert mode in PopulationDistribution.get_valid_modes()

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
		Get a dictionary mapping genome id to the path of their genome

		@param file_path_of_file_mapping_genome_id_to_paths: File path to file with format 'id \t path'
		@type file_path_of_file_mapping_genome_id_to_paths: str | unicode
		@param list_of_drawn_genome_id: List of genome identifiers
		@type list_of_drawn_genome_id: list[str|unicode]

		@return: genome ids mapped to their gnome file path
		@rtype: dict[str|unicode, str|unicode]
		"""
		mdt = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		mdt.read(file_path_of_file_mapping_genome_id_to_paths)
		genome_id_to_path_map = mdt.get_map(0, 1, unique_key=True)
		assert set(genome_id_to_path_map.keys()).issuperset(list_of_drawn_genome_id)
		return {genome_id: genome_id_to_path_map[genome_id] for genome_id in list_of_drawn_genome_id}

	def _move_genome_file(
		self, file_path_input, file_path_output, sequence_min_length, set_of_sequence_names=None, file_format="fasta"):
		"""
		Move genomes into project folder, cleaning it up in the process.
		Makes sure sequence ids are unique and descriptions/comments are removed
		Returns total length of all sequences

		@attention file_format: Anything but 'fasta' is not supported, yet

		@param file_path_input: File path to raw genome
		@type file_path_input: str | unicode
		@param file_path_output: Destination path
		@type file_path_output: str | unicode
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@param file_format: 'fasta' format by default.
		@type file_format: str | unicode

		@return: Total length of all sequences
		@rtype: int | long
		@raise Exception:
		"""
		assert self.validate_file(file_path_input)
		assert not self.validate_file(file_path_output, silent=True), "Overwriting files prohibited: '{}'".format(
			file_path_output)
		assert file_format == "fasta", "'{}' is not supported, yet.".format(file_format)
		if set_of_sequence_names is None:
			set_of_sequence_names = []
		with open(file_path_input, 'r') as stream_input, open(file_path_output, 'w') as stream_output:
			total_base_pairs = self._cleanup_and_filter_sequences(
				stream_input, stream_output, sequence_min_length, set_of_sequence_names, file_format)
		if total_base_pairs == 0:
			msg = "No valid sequences in '{}'".format(stream_input.name)
			self._logger.error(msg)
			raise Exception(msg)
		return total_base_pairs

	def _cleanup_and_filter_sequences(
		self, stream_input, stream_output, sequence_min_length, set_of_sequence_names, file_format="fasta"):
		"""

		@attention file_format: Anything but 'fasta' is not supported, yet

		@param stream_input: input stream of sequence file
		@type stream_input: file | FileIO | StringIO
		@param stream_output: Output stream
		@type stream_output: file | FileIO | StringIO
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@param file_format: 'fasta' format by default.
		@type file_format: str | unicode

		@return: Total length of all sequences (base pairs)
		@rtype: int | long
		"""
		total_base_pairs = 0
		for seq_record in SeqIO.parse(stream_input, file_format):
			# remove description, else art illumina messes up sam format
			seq_record.description = ''
			if len(seq_record.seq) < sequence_min_length:
				self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
					os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
				continue
			if seq_record.id in set_of_sequence_names:
				seq_record.id = self._get_new_name(seq_record.id, set_of_sequence_names)
			set_of_sequence_names.add(seq_record.id)
			# file_handler.write(">{}\n".format(sequence_id))
			# file_handler.writelines("{}\n".format(seq_record.seq))
			stream_output.write(seq_record.format(file_format))
			total_base_pairs += len(seq_record.seq)
		return total_base_pairs

	def move_genome_files(
		self, genome_id_to_path_map, directory_output, sequence_min_length, set_of_sequence_names=None):
		"""
		Move and clean up a list of genomes

		@param genome_id_to_path_map: Dictionary with file path by genome ids
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]
		@param directory_output: Output directory
		@type directory_output: str | unicode
		@param sequence_min_length: Minimum length of sequences
		@type sequence_min_length: int | long
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]

		@return:
		@rtype: dict[str|unicode, int|long]
		"""
		assert isinstance(genome_id_to_path_map, dict)
		if set_of_sequence_names is None:
			set_of_sequence_names = set()
		genome_id_to_total_length = {}
		for genome_id, genome_file_path in genome_id_to_path_map.iteritems():
			file_name = os.path.basename(genome_file_path)
			new_genome_file_path = os.path.join(directory_output, file_name)
			total_length = self._move_genome_file(
				genome_file_path, new_genome_file_path, sequence_min_length, set_of_sequence_names)
			genome_id_to_total_length[genome_id] = total_length
		return genome_id_to_total_length

	@staticmethod
	def _get_new_name(name, set_of_sequence_names):
		"""
		Get a unique sequence name

		@param name: Current sequence name
		@type name: str | unicode
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@return: Unique sequence name
		@rtype : str | unicode
		"""
		index = 0
		new_name = name
		while new_name in set_of_sequence_names:
			new_name = name + "_{}".format(index)
			index += 1
		return new_name

	def validate_format(self, list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
		"""
		Validate file format of a list of fasta files

		@param list_of_file_paths: List of fasta file paths
		@type list_of_file_paths: list[str|unicode]
		@param file_format: 'fasta' or 'fastq'
		@type file_format: str | unicode
		@param sequence_type: 'dna' or 'rna' or 'protein'
		@type sequence_type: str | unicode
		@param ambiguous: If true ambiguous characters are valid
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
	"""
		For the design of an artificial community
	"""
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
		@param column_name_genome_id: Column name of genome ids in the metadata table
		@type column_name_genome_id: str | unicode
		@param column_name_otu: Column name of otu ids in the metadata table
		@type column_name_otu: str | unicode
		@param column_name_novelty_category: Column name of novelty category in the metadata table
		@type column_name_novelty_category: str | unicode
		@param column_name_ncbi: Column name of taxonomic id assignment in the metadata table
		@type column_name_ncbi: str | unicode
		@param column_name_source: Column name of 'source' in the metadata table
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
		self, directory_output, genome_id_to_total_length, genome_id_to_abundance, genome_id_to_file_name,
		community_id):
		"""
		Write abundance file for each sample

		@param directory_output: Distribution files will be stored here
		@type directory_output: str | unicode
		@param genome_id_to_total_length: Total sequence length of each genome
		@type genome_id_to_total_length: dict[str|unicode, int|long]
		@param genome_id_to_abundance: Drawn distribution for each genome id
		@type genome_id_to_abundance: dict[str|unicode, list[float]]
		@param genome_id_to_file_name: genome file path of each genome id
		@type genome_id_to_file_name: dict[str|unicode, str|unicode]
		@param community_id: Id of a community
		@type community_id: int | str | unicode
		"""
		list_of_genome_id = genome_id_to_abundance.keys()
		distributions = len(genome_id_to_abundance[list_of_genome_id[0]])
		for index_sample in range(0, distributions):
			out_file_path = os.path.join(
				directory_output, self._filename_distribution_comunity.format(
					comunity_index=community_id, sample_index=index_sample))
			with open(out_file_path, 'w') as out_file_handler:
				for genome_id in genome_id_to_abundance:
					distribution = genome_id_to_abundance[genome_id][index_sample]
					file_name = genome_id_to_file_name[genome_id]
					total_length = genome_id_to_total_length[genome_id]
					out_file_handler.write("{id}\t{file_name}\t{dist}\t{len}\n".format(
						id=genome_id,
						file_name=file_name,
						dist=distribution,
						len=total_length))

	def design(
		self, community, number_of_samples, metadata_table,
		directory_out_distributions, directory_out_genomes, directory_out_metadata, directory_in_template=None,
		set_of_sequence_names=None, min_sequence_length=1):
		"""
		Design artificial communities, of a specific design, for a given number of samples

		@param community: Input data for the creation of a community
		@type community: Community
		@param number_of_samples: Amount of samples to be simulated
		@type number_of_samples: int
		@param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
		@type metadata_table: MetadataTable
		@param directory_out_distributions: Output directory where distribution files will be placed
		@type directory_out_distributions: str | unicode
		@param directory_out_genomes:  Output directory where cleaned up raw fasta files will be placed
		@type directory_out_genomes: str | unicode
		@param directory_out_metadata: Metadata tables of separated by chosen and not chosen genomes are written to here
		@type directory_out_metadata: str | unicode
		@param directory_in_template: contains template data for strain simulation
		@type directory_in_template: str | unicode
		@param set_of_sequence_names: Set of all previously used sequence names, making sure all will be unique
		@type set_of_sequence_names: set[str|unicode]
		@param min_sequence_length: Minimum sequence length
		@type min_sequence_length: int | long

		@return: List of drawn genome ids
		@rtype: list[str|unicode]
		"""
		assert isinstance(community, Community)
		assert isinstance(metadata_table, MetadataTable)

		if set_of_sequence_names is None:
			set_of_sequence_names = set()

		number_of_strains = community.genomes_total

		# pick how much a strain will be simulated
		genome_amounts = []
		strain_simulation = None
		if community.simulate_strains:
			strain_simulation = StrainSimulationWrapper(
				executable_sim=None,
				directory_template=directory_in_template,
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
		metadata_table_community = MetadataTable(logfile=self._logfile, verbose=self._verbose)
		metadata_table_community.read(community.file_path_metadata_table, column_names=True)
		strain_selector = StrainSelector(
			column_name_genome_id=self._column_name_genome_id,
			column_name_otu=self._column_name_otu,
			column_name_novelty_category=self._column_name_novelty_category,
			logfile=self._logfile, verbose=self._verbose, debug=self._debug
			)
		list_of_drawn_genome_id = strain_selector.get_drawn_genome_id(
			metadata_table=metadata_table_community,
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
		metadata_new_file_path = os.path.join(directory_out_metadata, new_file_name)
		metadata_table_community.write(
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
		metadata_table_community.reduce_rows_to_subset(list_of_drawn_genome_id, self._column_name_genome_id)
		metadata_table.concatenate(metadata_table_community, strict=False)

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
				meta_table=metadata_table,
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
			directory_output=directory_out_genomes,
			sequence_min_length=min_sequence_length,
			set_of_sequence_names=set_of_sequence_names)

		# write distribution file
		# genome_id_to_distributions = self._get_genome_id_to_distributions(list_of_drawn_genome_id, list_of_distributions)
		assert len(list_of_drawn_genome_id) == len(list_of_distributions)
		genome_id_to_distributions = dict(zip(list_of_drawn_genome_id, list_of_distributions))

		genome_id_to_file_name = self._get_genome_id_to_file_name(genome_id_to_path_map)
		self._write_distribution_files(
			directory_output=directory_out_distributions,
			genome_id_to_total_length=genome_id_to_total_length,
			genome_id_to_abundance=genome_id_to_distributions,
			genome_id_to_file_name=genome_id_to_file_name,
			community_id=community.id)
		return list_of_drawn_genome_id

	@staticmethod
	def _get_genome_id_to_file_name(genome_id_to_path_map):
		"""
		Extract file names for reuse in new folder

		@attention: Preserving original name might will cause issue if not unique

		@param genome_id_to_path_map: Dictionary with file path by genome ids
		@type genome_id_to_path_map: dict[str|unicode, str|unicode]

		@return: Map of genome id to a filename
		@rtype : dict[str|unicode, str|unicode]
		"""
		set_of_file_names = set()
		genome_id_to_file_name = {}
		for genome_id, file_path in genome_id_to_path_map.iteritems():
			filename = os.path.basename(file_path)
			assert filename not in set_of_file_names, "Filename '{}' is not unique!".format(set_of_file_names)
			set_of_file_names.add(filename)
			genome_id_to_file_name[genome_id] = filename
		return genome_id_to_file_name

	#####################
	#
	# merge communities
	#
	#####################

	def merge_communities(self, list_of_communities, directory_out_distributions, number_of_samples, metadata_table):
		"""
		Combine distributions of communities and adjust them according to their ratio.

		@param list_of_communities: List of community inputs
		@type list_of_communities: list[Community]
		@param directory_out_distributions: Output directory where distribution files will be placed
		@type directory_out_distributions: str | unicode
		@param number_of_samples: Amount of samples to be simulated
		@type number_of_samples: int | long
		@param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
		@type metadata_table: MetadataTable

		@return: List of combined distribution file paths for each sample
		@rtype: list[str|unicode]
		"""
		assert isinstance(list_of_communities, list)
		for community in list_of_communities:
			assert isinstance(community, Community)
		assert isinstance(metadata_table, MetadataTable)

		# read communities and adapt to ratio
		list_of_output_paths = []
		for index_sample in range(number_of_samples):
			communities = []
			list_of_community_total_abundance = [0] * len(list_of_communities)
			sample_total_abundance = 0

			# create list of community files for one sample
			list_of_file_paths = []
			for index_community, c in enumerate(list_of_communities):
				dist_file_name = os.path.join(
					directory_out_distributions,
					self._filename_distribution_comunity.format(
						comunity_index=index_community, sample_index=index_sample))
				list_of_file_paths.append(dist_file_name)

			metadata_table_community = MetadataTable(logfile=self._logfile, verbose=self._verbose)
			for index_community, file_path in enumerate(list_of_file_paths):
				communities.append(metadata_table_community.parse_file(file_path, column_names=False))
				# communities_length.append(0)

				genomes = set()
				for genome_id, file_name, abundance, genome_length in communities[index_community]:
					if genome_id in genomes:
						# print genome_id, filename, abundance, genome_length
						raise ValueError("Genome id '{}' not unique".format(genome_id))
					genomes.add(genome_id)
					list_of_community_total_abundance[index_community] += float(abundance)  # * float(sequence_info[4])
				sample_total_abundance += list_of_community_total_abundance[index_community]
				communities[index_community].close()

			# out.append(read_communities[0][0])
			list_of_community_factor = [0.0] * len(list_of_communities)
			for index_community, c in enumerate(list_of_file_paths):
				ratio = float(list_of_communities[index_community].ratio)
				community_total_abundance = float(list_of_community_total_abundance[index_community])
				current_proportion_in_sample = community_total_abundance / float(sample_total_abundance)
				list_of_community_factor[index_community] = ratio / current_proportion_in_sample
				# self.update_community(communities[index_community], factor)

			# join communities
			communities = []
			for index_community, file_path in enumerate(list_of_file_paths):
				communities.append(metadata_table_community.parse_file(file_path, column_names=False))

			# print_ratios(communities)
			file_path_output = os.path.join(
				directory_out_distributions, self._filename_distribution_comunity_joint.format(sample_index=index_sample))
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
		Stream out joined distribution for a sample

		@param communities: list of iterators over row cells from community abundance files
		@type communities: list[generator]
		@param list_of_community_factor: multiplication factor for each community to get right ratio
		@type list_of_community_factor: list[float]
		@param stream_output: joined distribution information output
		@type stream_output: file | FileIO | StringIO | basestring
		"""
		line_format = "{gid}\t{filename}\t{distr}\t{length}\n"
		for community_index, community in enumerate(communities):
			factor = list_of_community_factor[community_index]
			for genome_id, filename, abundance, genome_length in community:
				stream_output.write(line_format.format(
					gid=genome_id,
					filename=filename,
					distr=float(abundance) * factor,
					length=genome_length,
				))
			community.close()
