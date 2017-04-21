""" Collection of handy tools used by various modules of the CONSTANd++ workflow. """


def unnest(x):
	""" returns un-nested version of level 1 nested list x."""
	return [e for sublist in x for e in sublist]
