#!/usr/bin/python

__author__ = 'hofmann'

import unittest
import os
import sys
import shutil
from communitydesign import CommunityDesign, Community
from scripts.GenomePreparation.genomepreparation import GenomePreparation
from taxonomicprofile import TaxonomicProfile
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy

# ##################################
#
#          PrepareStrains
#
# ##################################


class DefaultSetUpPrepareStrains(unittest.TestCase):
	_test_case_id = 0
	_success = False
	log_filename = 'unittest_log.txt'
	filename_metadata = "meta_data.tsv"

	dir_input = "unittest_input"
	dir_output = "unittest_output_PS_{}"

	filename_input = "path_genoms.tsv"

	def setUp(self):
		self.dir_output = self.dir_output.format(self._test_case_id)
		if os.path.isdir(self.dir_output):
			shutil.rmtree(self.dir_output)
		os.mkdir(self.dir_output)
		sys.stderr.write("\n{}... ".format(self._test_case_id))
		DefaultSetUpPrepareStrains._test_case_id += 1

		logfile = os.path.join(self.dir_output, self.log_filename)
		self.file_stream = open(logfile, 'a')
		self.test_object = GenomePreparation(
			logfile=self.file_stream,
			verbose=False,
			debug=True)

		self.number_of_strains_per_otu = 3

	def tearDown(self):
		self.test_object = None
		self.file_stream.close()
		self.file_stream = None
		if self._success:
			shutil.rmtree(self.dir_output)


class TestMethodsPrepareStrains(DefaultSetUpPrepareStrains):

	def test_move_genome_files(self):
		list_of_sequence_names = set()
		list_of_drawn_genome_id = ["126738.1", "1071768.1", "1244532.1", "234267.1", "443255.1"]
		min_length = 1
		directory_genomes = self.dir_output
		file_path_genome_id_to_path_map = os.path.join(self.dir_input, self.filename_input)
		genome_id_to_path_map = self.test_object._get_genome_id_to_path_map(
			file_path_genome_id_to_path_map, list_of_drawn_genome_id)
		dict_genome_id_to_path = self.test_object.move_genome_files(
			genome_id_to_path_map, directory_genomes, min_length, list_of_sequence_names)
		self.assertIsInstance(dict_genome_id_to_path, dict, "{}".format(dict_genome_id_to_path))
		self._success = True


# ##################################
#
#          CommunityDesign
#
# ##################################


class DefaultSetUpCommunityDesign(unittest.TestCase):
	_test_case_id = 0
	_success = False
	log_filename = 'unittest_log.txt'
	filename_metadata = "meta_data.tsv"
	filename_metadata2 = "meta_data2.tsv"

	dir_input = "unittest_input"
	dir_output = "unittest_output_CD_{}"

	def setUp(self):
		self.dir_output = self.dir_output.format(self._test_case_id)
		if os.path.isdir(self.dir_output):
			shutil.rmtree(self.dir_output)
		os.mkdir(self.dir_output)
		sys.stderr.write("\n{}... ".format(self._test_case_id)),
		DefaultSetUpCommunityDesign._test_case_id += 1

		logfile = os.path.join(self.dir_output, self.log_filename)
		self.file_stream = open(logfile, 'a')
		self.test_object = CommunityDesign(
			max_processors=4,
			tmp_dir=self.dir_output,
			logfile=self.file_stream,
			verbose=False,
			debug=False,
			seed="seed")

		self.taxonomy_directory = os.path.join(self.dir_input, "ncbi-taxonomy_20150130")
		self.filename_path_genoms = "path_genoms.tsv"
		self.filename_path_gff = "path_genoms_gen_annotation.tsv"
		self.directory_distributions = os.path.join(self.dir_output, "distr")
		self.directory_genomes = os.path.join(self.dir_output, "source")
		self.directory_metadata = os.path.join(self.dir_output, "meta")
		os.mkdir(self.directory_distributions)
		os.mkdir(self.directory_genomes)
		os.mkdir(self.directory_metadata)

	def tearDown(self):
		self.test_object = None
		self.file_stream.close()
		self.file_stream = None
		if self._success:
			shutil.rmtree(self.dir_output)


class TestMethodsCommunityDesign(DefaultSetUpCommunityDesign):

	def test_get_drawn_genome_id(self):
		expected_result = [
			'simulated_443255.1.Taxon008',
			'simulated_1244532.1.Taxon029',
			'234267.1',
			'simulated_1244532.1.Taxon036',
			'simulated_234267.1.Taxon014',
			'1244532.1',
			'simulated_234267.1.Taxon005',
			'simulated_234267.1.Taxon007',
			'443255.1',
			'1071768.1']
		# file_path_metadata_table = os.path.join(self.dir_input, self.filename_metadata)
		metadata_table_comby = MetadataTable(verbose=False)
		file_path_metadata_table = os.path.join(self.dir_input, self.filename_metadata)
		file_path_genome_locations = os.path.join(self.dir_input, self.filename_path_genoms)
		file_path_gff_locations = os.path.join(self.dir_input, self.filename_path_gff)

		community = Community(
			identifier=str(1),
			genomes_total=10,
			genomes_real=4,
			limit_per_otu=2,
			file_path_metadata_table=file_path_metadata_table,
			file_path_genome_locations=file_path_genome_locations,
			file_path_gff_locations=file_path_gff_locations,
			ratio=1,
			mode="differential",
			log_mu=1,
			log_sigma=2,
			gauss_mu=1,
			gauss_sigma=3,
			verbose=False)

		set_of_sequence_names = set()
		list_of_drawn_genome_id = self.test_object.design(
			community,
			number_of_samples=5,
			metadata_table=metadata_table_comby,
			directory_out_distributions=self.directory_distributions,
			directory_out_genomes=self.directory_genomes,
			directory_out_metadata=self.directory_metadata,
			directory_in_template=None,
			set_of_sequence_names=set_of_sequence_names,
			min_sequence_length=1)
		# print list_of_drawn_genome_id
		self.assertListEqual(list_of_drawn_genome_id, expected_result)
		self._success = True

	def test_merge_communities(self):
		metadata_table_comby = MetadataTable(verbose=False)
		file_path_metadata_table = os.path.join(self.dir_input, self.filename_metadata)
		file_path_metadata_table2 = os.path.join(self.dir_input, self.filename_metadata2)
		file_path_genome_locations = os.path.join(self.dir_input, self.filename_path_genoms)
		file_path_gff_locations = os.path.join(self.dir_input, self.filename_path_gff)

		community0 = Community(
			identifier=str(0),
			genomes_total=10,
			genomes_real=4,
			limit_per_otu=2,
			file_path_metadata_table=file_path_metadata_table,
			file_path_genome_locations=file_path_genome_locations,
			file_path_gff_locations=file_path_gff_locations,
			ratio=1,
			mode="differential",
			log_mu=1,
			log_sigma=2,
			gauss_mu=1,
			gauss_sigma=3,
			verbose=False)
		community1 = Community(
			identifier=str(1),
			genomes_total=10,
			genomes_real=1,
			limit_per_otu=2,
			file_path_metadata_table=file_path_metadata_table2,
			file_path_genome_locations=file_path_genome_locations,
			file_path_gff_locations=file_path_gff_locations,
			ratio=1,
			mode="differential",
			log_mu=1,
			log_sigma=2,
			gauss_mu=1,
			gauss_sigma=3,
			verbose=False)

		number_of_samples = 5
		list_of_sequence_names = set()
		list_of_drawn_genome_id = self.test_object.design(
			community0,
			number_of_samples=5,
			metadata_table=metadata_table_comby,
			directory_out_distributions=self.directory_distributions,
			directory_out_genomes=self.directory_genomes,
			directory_out_metadata=self.directory_metadata,
			directory_in_template=None,
			set_of_sequence_names=list_of_sequence_names,
			min_sequence_length=1)

		list_of_drawn_genome_id.extend(self.test_object.design(
			community1,
			number_of_samples=5,
			metadata_table=metadata_table_comby,
			directory_out_distributions=self.directory_distributions,
			directory_out_genomes=self.directory_genomes,
			directory_out_metadata=self.directory_metadata,
			directory_in_template=None,
			set_of_sequence_names=list_of_sequence_names,
			min_sequence_length=1))

		list_of_communities = [community0, community1]
		list_of_output_paths = self.test_object.merge_communities(
			list_of_communities,
			directory_out_distributions=self.directory_distributions,
			number_of_samples=number_of_samples,
			metadata_table=metadata_table_comby)
		# print list_of_drawn_genome_id
		self.assertEqual(len(list_of_output_paths), number_of_samples)
		taxonomy = NcbiTaxonomy(taxonomy_directory=self.taxonomy_directory, verbose=False)
		tp = TaxonomicProfile(taxonomy, logfile=None, verbose=False, debug=False)
		tp.write_taxonomic_profile_from_abundance_files(
			metadata_table=metadata_table_comby,
			list_of_file_paths=list_of_output_paths,
			directory_output=self.dir_output
		)
		out_metadata = os.path.join(self.dir_output, self.filename_metadata)
		metadata_table_comby.write(out_metadata)
		self._success = True


if __name__ == '__main__':
	suite0 = unittest.TestLoader().loadTestsFromTestCase(TestMethodsPrepareStrains)
	suite1 = unittest.TestLoader().loadTestsFromTestCase(TestMethodsCommunityDesign)
	all_tests = unittest.TestSuite([suite0, suite1])
	unittest.TextTestRunner(verbosity=2).run(all_tests)
	# unittest.TextTestRunner(verbosity=2).run(suite0)
