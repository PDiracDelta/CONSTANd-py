from unittest import TestCase
from dataIO import *
from constand import constand
from numpy.testing import assert_array_almost_equal


class TestConstand(TestCase):
	def setUp(self):
		self.file_in = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/data/MB_noapostrophes.tsv'  # TEST
		self.delim_in = '\t'
		self.header_in = None
		self.accuracy = 1e-2
		self.maxIterations = 50
		self.path_out = '/home/pdiracdelta/Documents/KUL/Master of Bioinformatics/Thesis/data/MB_result.tsv'  # TEST
		self.delim_out = '\t'
		self.intensities, self.df = importDataFrame(self.file_in, self.delim_in)

	def test_constand_result(self):
		correctResult = np.asarray(importDataFrame(file_in=self.path_out, header=self.header_in))
		normalizedIntensities, convergenceTrail, R, S = constand(self.intensities, self.accuracy, self.maxIterations)
		assert_array_almost_equal(normalizedIntensities, correctResult,
		                          err_msg="CONSTANd result inconsistent with earlier calculations.")
		assert_array_almost_equal(((self.intensities * S).T * R).T, correctResult,
		                          err_msg="R and S matrices do not produce correct R*A*S result.")
		self.assertTrue(np.all(convergenceTrail[0:-2] > self.accuracy) and convergenceTrail[-1] < self.accuracy,
		                "The residual errors in the convergenceTrail are inconsistent with the imposed accuracy.")
