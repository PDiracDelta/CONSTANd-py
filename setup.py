from setuptools import setup, find_packages
from qcquan import __version__ as version

setup(name='qcquan',
      version=version,
      description='QCQuan workflow',
      url='https://constand.vito.be',
      author='Joris Van Houtven',
      author_email='joris.vanhoutven@uhasselt.be',
      license='VITO',
      packages=find_packages(),
      install_requires=['werkzeug', 'flask', 'requests', 'scipy', 'pandas', 'numpy', 'sklearn',
                        'statsmodels', 'weasyprint', 'matplotlib', 'wtforms'],
      zip_safe=False)