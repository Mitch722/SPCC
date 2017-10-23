

class Data(object):
	"""docstring for Data"""
	def __init__(self, M, disB='g', no_subset):
		super(Data, self).__init__()
		self.M = M
		self.disB = disB
		self.no_subset = no_subset
		self.n = self.M / self.no_subSet


