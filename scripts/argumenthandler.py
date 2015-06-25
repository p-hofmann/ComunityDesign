__author__ = 'hofmann'

import os
import sys
import argparse
import tempfile
from scripts.configparserwrapper import ConfigParserWrapper
from scripts.Validator.validator import Validator


class ArgumentHandler(Validator):
	"""Reading pipeline configuration from file and from passed arguments"""

	_label = "ArgumentHandler"

	_seed = None  # "00000000000000000000000000000000000000000000000000"
	_pooled_gsa = False

	_phase = 0
	_debug = False
	_verbose = False

	# ############
	# [main]
	# ############
	_tmp_dir = None
	_directory_ncbi_taxdump = None
	_directory_pipeline = None
	_file_path_config = None
	_directory_output = None
	_max_processors = 1
	_dataset_id = ''
	# TODO: read from config!
	_file_path_samtools = None  # "samtools"

	# ############
	# [read_simulator]
	# ############
	_sample_size_in_base_pairs = None
	_base_pairs_multiplication_factor = float(1000000000)  # 10**9

	_read_simulator_type = None
	_error_profile = None
	_fragment_size_standard_deviation_in_bp = None
	_fragments_size_mean_in_bp = None

	# ############
	# [sampledesign]
	# ############
	_strain_simulation_template = "tools/sgEvolver/simulation_dir"
	_number_of_samples = None
	_file_path_plasmid_sequence_names = None

	# ############
	# [community*]
	# ############
	_number_of_communities = None
	_communities = []

	# ############
	# [comdesign]
	# ############
	_abundance_sigma = None  # 2
	_abundance_mean = None  # 1
	# _use_strain_simulation = None
	_abundance_file_prefix = "distr"

	# ############
	# [differential]
	# ############

	# multiple_samples_modus = None
	_view_distribution = False  # better [comdesign] ??

	# subfolder_names
	_folder_name_comunity_design = "comunity_design"
	_folder_name_distribution = "distributions"
	_folder_name_genomes = "source_genomes"
	_folder_name_meta_data = "meta_data"
	# folder_name_simulated = "simulated_genomes"
	_folder_name_bam = "bam"
	_folder_name_sam = "sam"
	_folder_name_fastq = "fastq"
	_folder_name_logfiles = "logfiles"
	_sub_folders_sample = [_folder_name_bam, _folder_name_sam, _folder_name_fastq, _folder_name_logfiles]
	_sub_folders_output = [_folder_name_distribution, _folder_name_genomes]

	# filenames
	_ncbi_ref_files = ["nodes.dmp", "merged.dmp", "names.dmp"]

	_filename_anonymous_reads = "anonymous_reads.fq"
	# filename_reads_anonymous_mapping = "reads_anonymous_mapping.tsv"

	_filename_anonymous_gsa = "anonymous_gsa.fasta"
	# filename_gsa_anonymous_mapping = "gsa_anonymous_mapping.tsv"

	_filename_anonymous_gsa_pooled = "anonymous_gsa_pooled.fasta"
	# filename_pooled_gsa_mapping = "pooled_" + filename_gsa_anonymous_mapping

	# filename_gsa = "gsa.fasta"
	# filename_pooled_gsa = "pooled_" + filename_gsa
	_filename_reads_mapping = "reads_mapping.tsv"
	_filename_gsa_mapping = "gsa_mapping.tsv"
	_filename_pooled_gsa_mapping = "gsa_pooled_mapping.tsv"

	_filename_log = "pipeline.log"
	_filename_metadata = "meta_data.tsv"

	def __init__(self, args=None):
		self._directory_pipeline = self._get_directory_pipeline()

		options = self._get_parser_options(args)
		super(ArgumentHandler, self).__init__(options.logfile, options.verbose)

		self._valid_arguments = True
		self._read_options(options)
		if not self._valid_arguments:
			return
		self._config = ConfigParserWrapper(self._file_path_config)
		self._read_config()
		if not self._valid_arguments:
			return
		self._check_values()
		if not self._valid_arguments:
			return
		if self._directory_ncbi_taxdump is not None and self._phase < 2:
			assert self.validate_dir(self._directory_ncbi_taxdump, file_names=self._ncbi_ref_files)
		options = ArgumentHandler(args)
		if options.tmp_dir is None:
			options.tmp_dir = tempfile.gettempdir()
		assert self.validate_dir(options.tmp_dir)

	def _get_directory_pipeline(self):
		return self.get_full_path(os.path.dirname(os.path.realpath(sys.argv[0])))

	def get_logfile(self):
		return self._logfile

	def to_file(self, file_path):
		file_directory = os.path.dirname(file_path)
		if not os.path.isdir(file_directory):
			self._logger.error("Directory does not exist: '{}'".format(file_directory))
			return
		with open(file_path, 'w') as file_handler:
			file_handler.write(self.to_string())

	def to_string(self):
		expected_output_size = self._expected_output_size_in_giga_byte()
		if expected_output_size < 0.001:
			expected_output_size = 0.001
		stages = ["Full run through", "Community creation only", "Read simulator only"]
		result_string = """
[Parameter]:
Dataset id:\t\t{dataset}

_main_
Max Processes:\t\t{pool}
Pipeline directory:\t'{pipe}'
Config file:\t\t'{config}'
Output directory:\t'{out}'
Phase:\t\t\t{stage}
Expected output size:\t< {out_size:.3f} GigaByte

_sample_design_
Number of samples:\t{samples}
Pooled GSA:\t\t{pgsa}

_read_simulator_
Expected total size:\t{bps} base pairs / {gbps:.2f} gbp
Read simulator:\t\t'{readsim}'
Error profile:\t\t'{error}'
Fragment size:
	Mean:\t\t{fmean}
	Standard deviation:\t{fsd}


_comunity_design_
""".format(
			config=self._file_path_config,
			pipe=self._directory_pipeline,
			out=self._directory_output,
			stage=stages[self._phase],
			out_size=expected_output_size,
			bps=self._sample_size_in_base_pairs,
			gbps=float(self._sample_size_in_base_pairs)/self._base_pairs_multiplication_factor,
			samples=self._number_of_samples,
			readsim=self.read_simulator,
			error=self._error_profile,
			fsd=self._fragment_size_standard_deviation_in_bp,
			fmean=self._fragments_size_mean_in_bp,
			dataset=self._dataset_id,
			pool=self._max_processors,
			pgsa=str(self._pooled_gsa)
			# plasmid=self.plasmid_file
			)

		for index, community in enumerate(self._communities):
			com_string = """
community{i}:
	number of genomes = {ngs}
	metadata file = {mf}
	genome id to file = {gtf}
	simulated ration = {evolve}
	number per otu = {npo}
	ratio = {ratio}
	modus = {modus}
	mu = {mu}
	sigma = {sigma}
	view = {view}
""".format(
				i=index,
				ngs=community['num_genomes'],
				evolve=community['evolve'],
				npo=community['num_per_otu'],
				mf=community['metadata_file'],
				gtf=community['genomes_to_file'],
				ratio=community['ratio'],
				modus=community['modus'],
				# sn=community['sample_number'],
				mu=community['mu'],
				sigma=community['sigma'],
				view=community['view']
			)

			result_string += com_string
		return result_string

	def is_valid(self):
		return self._valid_arguments

	def _check_values(self):
		if self._max_processors is None:
			self._max_processors = 1

		if self._tmp_dir is None:
			self._tmp_dir = tempfile.gettempdir()

		subfolders = ["scripts", "tools"]

		if self._directory_pipeline is None:
			self._logger.error("Pipeline directory is required!")
			self._valid_arguments = False
			return
		elif not self.validate_dir(
			self._directory_pipeline,
			sub_directories=subfolders, key="pipeline directory"):
			self._valid_arguments = False
			return

		if self._directory_output is None:
			self._logger.error("'-o' Output directory is required!")
			self._valid_arguments = False
			return

		self._directory_output = self.get_full_path(self._directory_output)

		if not self.validate_dir(self._directory_output, only_parent=True, key='-o'):
			self._valid_arguments = False
			return

		if self._sample_size_in_base_pairs is None:
			self._logger.error("'-bp' A size in basepairs must be given!")
			self._valid_arguments = False
			return
		elif not self.validate_number(self._sample_size_in_base_pairs, minimum=1, key='-bp'):
			self._valid_arguments = False
			return

		if self._number_of_samples is None:
			self._logger.error("'-ns' No number of samples given!")
			self._valid_arguments = False
			return

		if not self.validate_number(self._number_of_samples, minimum=1, key='-ns'):
			self._valid_arguments = False
			return

		for index in range(self._number_of_communities):
			file_list = [self._communities[index]['metadata_file'], self._communities[index]['genomes_to_file']]
			for file_path in file_list:
				if file_path is None:
					self._logger.error("[community{index}] File missing!".format(index=index))
					self._valid_arguments = False
					return
				key = "community{index}".format(index=index)
				if not self.validate_file(file_path, key=key):
					self._valid_arguments = False
					return

		if not self.validate_number(self._number_of_samples, maximum=10):
			if not self._verbose:
				self._logger.error("{} samples are an unusual amount. User input required!".format(self._number_of_samples))
				self._valid_arguments = False
				return

			message = "{} samples are an unusual amount. Are you sure you want to continue? [y/n]".format(
				self._number_of_samples)
			if not self.get_confirmation(message):
				self._valid_arguments = False
				return

		if self.read_simulator is None:
			self._logger.error("'-rs' No read simulator given!")
			self._valid_arguments = False
			return

		if self.read_simulator == 'art':
			if self._error_profile is None:
				self._logger.error("'-ep' An error profile for 'art' is needed: 'mi' or 'hi'!")
				self._valid_arguments = False
				return

			if self._fragments_size_mean_in_bp is None:
				self._logger.error("'-fmean' For the simulation with 'art' a mean size of the fragments is required!")
				self._valid_arguments = False
				return
			elif not self.validate_number(self._fragments_size_mean_in_bp, minimum=1, key='-fmean'):
				self._valid_arguments = False
				return

			if self._fragment_size_standard_deviation_in_bp is None:
				self._logger.error("'-fsd' For the simulation with 'art' a standard_deviation of the fragments size is required!")
				self._valid_arguments = False
				return
			elif not self.validate_number(self._fragment_size_standard_deviation_in_bp, minimum=1, key='-fsd'):
				self._logger.error(
					"'-fsd' The standard_deviation of the fragments size must be a positive number: '{}'".format(
						self._fragment_size_standard_deviation_in_bp))
				self._valid_arguments = False
				return

		if self._file_path_plasmid_sequence_names is not None and not self.validate_file(
			self._file_path_plasmid_sequence_names, key=''):
			self._valid_arguments = False
			return

		expected_output_size = self._expected_output_size_in_giga_byte()
		expected_tmp_size = expected_output_size / self._number_of_samples
		directory_out = self._directory_output
		directory_tmp = self._tmp_dir
		if not os.path.isdir(directory_out):
			directory_out = os.path.dirname(directory_out)

		user_input_required = False
		if not self.validate_free_space(directory_tmp, required_space_in_gb=expected_tmp_size):
			user_input_required = True
		elif not self.validate_free_space(directory_out, required_space_in_gb=expected_output_size):
			user_input_required = True
		elif expected_output_size > 100:
			# message = "The output will require approximately {} GigaByte.".format(expected_output_size)
			self._logger.warning("The output will require approximately {} GigaByte.".format(expected_output_size))
		if user_input_required:
			if not self._verbose:
				self._logger.error("Continuation only possible with enabled user input!>")
				self._valid_arguments = False
				return

			message = "Are you sure you want to continue? [y/n]"
			if not self.get_confirmation(message):
				self._valid_arguments = False
				return

		if self._dataset_id is None:
			self._dataset_id = ''

		if self._phase is None:
			self._phase = 0

		if self._phase > 1 and not self.validate_dir(self._directory_output, key='-o'):
			self._valid_arguments = False
			return

		if self._phase == 2:
			sub_folders = [self._folder_name_distribution, self._folder_name_genomes]
			if not self.validate_dir(self._directory_output, sub_directories=sub_folders, key='-o'):
				self._valid_arguments = False
				return

		if self._phase < 2:
			if self._directory_ncbi_taxdump is None:
				self._logger.error("NCBI taxdump directory is required!")
				self._valid_arguments = False
				return

			self._directory_ncbi_taxdump = self.get_full_path(self._directory_ncbi_taxdump)
			if not self.validate_dir(self._directory_ncbi_taxdump, file_names=self._ncbi_ref_files, key='-ncbi'):
				self._valid_arguments = False
				return
			return

		# TODO: check other values
		return

	def _read_config(self):
		sections = ['main', 'readsimulator', 'sampledesign']
		invalid_sections = self._config.validate_sections(sections)
		if invalid_sections is not None:
			self._logger.error("Missing section '{}' in the configuration file.".format(", ".join(invalid_sections)))
			self._valid_arguments = False
			return

		# ##########
		# [main]
		# ##########

		if self._phase is None:
			self._phase = self._config.get_value("main", "phase", is_digit=True)

		if self._max_processors is None:
			self._max_processors = self._config.get_value("main", "max_processes", is_digit=True)

		if self._dataset_id is None:
			self._dataset_id = self._config.get_value("main", "dataset_id")

		if self._dataset_id is None:
			self._dataset_id = self._config.get_value("main", "dataset_id", is_digit=True)

		if self._directory_output is None:
			self._directory_output = self._config.get_value("main", "output_directory")

		if self._tmp_dir is None:
			self._tmp_dir = self._config.get_value("main", "temp_directory")

		if self._sample_size_in_base_pairs is None:
			config_value = self._config.get_value("sampledesign", "size", is_digit=True)
			if config_value is not None:
				self._sample_size_in_base_pairs = long(config_value * self._base_pairs_multiplication_factor)

		# ##########
		# [sampledesign]
		# ##########

		if self._directory_ncbi_taxdump is None:
			self._directory_ncbi_taxdump = self._config.get_value("sampledesign", "ncbi_taxdump_directory")

		if self._pooled_gsa is None:
			self._pooled_gsa = self._config.get_value("sampledesign", "pooled_gsa", is_boolean=True)

		if self._number_of_samples is None:
			self._number_of_samples = self._config.get_value("sampledesign", "number_of_samples", is_digit=True)

		if self._number_of_communities is None:
			self._number_of_communities = self._config.get_value('sampledesign', 'num_communities', is_digit=True)
		if self._number_of_communities is None:
			self._logger.error("Bad number of communities!")
			self._valid_arguments = False
			return

		for index_community in range(self._number_of_communities):
			if self._config.validate_sections(["community{}".format(index_community)]) is not None:
				self._logger.error("Missing 'community{}' section in config file".format(index_community))
				self._valid_arguments = False
				return
			community_name = "community{}".format(index_community)
			self._communities.append({})
			self._communities[index_community]['num_genomes'] = self._config.get_value(
				community_name, 'num_genomes', is_digit=True)
			self._communities[index_community]['num_real_genomes'] = self._config.get_value(
				community_name, 'num_real_genomes', is_digit=True)
			self._communities[index_community]['num_per_otu'] = self._config.get_value(
				community_name, 'number_per_otu', is_digit=True)
			self._communities[index_community]['metadata_file'] = self._config.get_value(community_name, 'metafile')
			self._communities[index_community]['genomes_to_file'] = self._config.get_value(community_name, 'genomes_to_file')
			self._communities[index_community]['ratio'] = self._config.get_value(community_name, 'ratio', is_digit=True)
			self._communities[index_community]['mu'] = self._config.get_value(community_name, 'mu', is_digit=True)
			self._communities[index_community]['sigma'] = self._config.get_value(community_name, 'sigma', is_digit=True)
			# self.communities[i]['sample_number'] = self._config.get_value(community_name, 'sample_number', True)
			self._communities[index_community]['view'] = self._config.get_value(community_name, 'view', is_boolean=True)
			self._communities[index_community]['modus'] = self._config.get_value(community_name, 'modus')
			evolve = self._config.get_value(community_name, 'evolved_rate', is_digit=True)
			if evolve is None:
				evolve = 0
			self._communities[index_community]['evolve'] = evolve

		# ##########
		# [sampledesign]
		# ##########

		if self.read_simulator is None:
			self.read_simulator = self._config.get_value("readsimulator", "type")

		if self._error_profile is None:
			self._error_profile = self._config.get_value("readsimulator", "profile")

		if self._fragment_size_standard_deviation_in_bp is None:
			self._fragment_size_standard_deviation_in_bp = self._config.get_value(
				"readsimulator", "fragment_size_standard_deviation", is_digit=True)

		if self._fragments_size_mean_in_bp is None:
			self._fragments_size_mean_in_bp = self._config.get_value("readsimulator", "fragments_size_mean", is_digit=True)

	def _expected_output_size_in_giga_byte(self):
		# very rough estimation
		expected_output_size = self._sample_size_in_base_pairs / float(1000000000) * 3 * 2
		expected_output_size *= self._number_of_samples
		return expected_output_size

	def _read_options(self, options):
		if not self.validate_file(options.config_file):
			self._valid_arguments = False
			return
		self._file_path_config = self.get_full_path(options.config_file)
		self._phase = options.phase
		self._verbose = options.verbose
		self._debug = options.debug_mode
		self._dataset_id = options.data_set_id
		self._max_processors = options.max_processors
		self._directory_output = options.output_directory
		self._sample_size_in_base_pairs = options.sample_size_gbp
		if self._sample_size_in_base_pairs is not None:
			self._sample_size_in_base_pairs = long(options.sample_size_gbp * self._base_pairs_multiplication_factor)
		self.read_simulator = options.read_simulator
		self._error_profile = options.error_profile
		self._fragment_size_standard_deviation_in_bp = options.fragment_size_standard_deviation
		self._fragments_size_mean_in_bp = options.fragments_size_mean
		# self.plasmid_file = options.plasmid_file
		self._number_of_samples = options.number_of_samples
		self._pooled_gsa = options.pooled_gsa

	@staticmethod
	def _get_parser_options(args=None):
		parser = argparse.ArgumentParser(
			description="Pipeline for the simulation of a metagenome",
			formatter_class=argparse.RawTextHelpFormatter)
		parser.add_argument(
			"-id", "--data_set_id",
			default=None,
			help="Id of the dataset, part of prefix of read/contig sequence ids!")
		parser.add_argument(
			"-p", "--max_processors",
			default=None,
			type=int,
			help="Maximum number processors used in parallel")
		parser.add_argument(
			"-verbose", "--verbose",
			action='store_true',
			default=None,
			help="Display more information!")
		parser.add_argument(
			"-debug", "--debug_mode",
			action='store_true',
			default=None,
			help="Activate debug modus. Temporary folders will not be deleted!")
		parser.add_argument(
			"-s", "--phase",
			default=None,
			type=int,
			choices=[0, 1, 2],
			help='''Available options: 0,1,2. Default: 0
0 -> Full run through,
1 -> Only Comunity creation,
2 -> Only Readsimulator
''')
		parser.add_argument(
			"-log", "--logfile",
			default=None,
			type=str,
			help="File path. Output will also written to this log file")

		# ##########
		# i/o
		# ##########

		group_input = parser.add_argument_group('i/o')
		group_input.add_argument(
			"config_file",
			type=str,
			help="Path to the configuration file of the pipeline")
		group_input.add_argument(
			"-o", "--output_directory",
			type=str,
			default=None,
			help="Directory in which the simulated data will be put into")

		# ##########
		# read simulator
		# ##########

		group_read_simulator = parser.add_argument_group('read simulator')
		group_read_simulator.add_argument(
			"-rs", "--read_simulator",
			type=str,
			default=None,
			choices=['art', 'pirs', 'pacbio'],
			help='''choose read simulator:
	'pacbio',
	'art' (Illumina),
	'pirs' (Illumina).
Currently only 'art' is supported.''')
		group_read_simulator.add_argument(
			"-ep", "--error_profile",
			type=str,
			default=None,
			choices=['mi', 'hi', 'hi150'],
			help='''Art Illumina error profiles:
	'mi': MiSeq 250bp,
	'hi': HiSeq 100bp,
	'hi150': HiSeq 150bp''')
		group_read_simulator.add_argument(
			"-bp", "--sample_size_gbp",
			type=float,
			default=None,
			help="The approximate total size the output will have in giga base pairs")
		group_read_simulator.add_argument(
			"-fmean", "--fragments_size_mean",
			type=int,
			default=None,
			help='''Mean size of the fragments in base pairs.''')
		group_read_simulator.add_argument(
			"-fsd", "--fragment_size_standard_deviation",
			type=int,
			default=None,
			help="Standard deviation from the mean size of fragments in base pairs.")

		# ##########
		# sample design
		# ##########

		group_community_design = parser.add_argument_group('sample design')
		group_community_design.add_argument(
			"-ns", "--number_of_samples",
			type=int,
			default=None,
			help='''Number of samples to be made''')
		group_community_design.add_argument(
			"-pgsa", "--pooled_gsa",
			default=None,
			action='store_true',
			help="Generate a pooled gsa from all samples")

		if args is None:
			return parser.parse_args()
		else:
			return parser.parse_args(args)

	def get_confirmation(self, message):
		user_input = raw_input("{}\n>".format(message)).lower()
		while True:
			if self.is_boolean_state(user_input):
				return self.get_boolean_state(user_input)
			user_input = raw_input("Please type 'n' for no, or 'y' for yes:\n>").lower()
