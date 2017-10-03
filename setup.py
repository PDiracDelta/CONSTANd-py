from setuptools import setup, find_packages

setup(name='qcquan',
      version='0.2',
      description='QCQuan workflow',
      url='https://constand.vito.be',
      author='Joris Van Houtven',
      author_email='joris.vanhoutven@uhasselt.be',
      license='VITO',
      packages=find_packages(),
      install_requires=['werkzeug', 'flask', 'requests', 'scipy', 'pandas', 'numpy', 'sklearn',
                        'statsmodels', 'weasyprint', 'matplotlib', 'wtforms'],
      zip_safe=False)