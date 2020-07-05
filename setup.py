from setuptools import setup, find_packages

from constand.constand import __version__

setup(name='CONSTANd',
      version=__version__,
      description='CONSTANd normalization method for relative quantification of data matrices based on the RAS/IPFP procedure.',
      url='https://qcquan.net/constand',
      author='Joris Van Houtven',
      author_email='joris.vanhoutven@uhasselt.be',
      license='VITO',
      packages=find_packages(),
      install_requires=['numpy'],
      zip_safe=False)