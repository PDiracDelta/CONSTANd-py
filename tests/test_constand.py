from unittest import TestCase
from dataIO import *
from constand import constand


class TestConstand(TestCase):
	def setUp(self):
		self.path_in = '../data/MB_noapostrophes.tsv'  # TEST
		self.delim_in = '\t'
		self.accuracy = 1e-2
		self.maxIterations = 50
		self.path_out = '../data/MB_result.tsv'  # TEST
		self.delim_out = '\t'
		self.intensities, self.df = importData(self.path_in, self.delim_in)

	def test_constand_result(self):
		correctResult, = importData(self.path_out)
		normalizedIntensities, convergenceTrail, R, S = constand(self.intensities, self.accuracy, self.maxIterations)
		self.assertEqual(normalizedIntensities, correctResult,
		                 "CONSTANd result inconsistent with earlier calculations.")
		self.assertEqual(R * self.intensities * S, correctResult,
		                 "R and S matrices do not produce correct R*A*S result.")
		self.assertTrue(convergenceTrail[0:-2] > self.accuracy > convergenceTrail[-1],
		                "The residual errors in the convergenceTrail are inconsistent with the imposed accuracy.")
