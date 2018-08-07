from setuptools import setup
from setuptools import find_packages

package_name = 'scikit-discovery'

package_list = find_packages()

with open("README.md", 'r', encoding='utf-8') as rfile:
    readme = rfile.read()

setup(name     = package_name,
      version  = '0.9.18',
      packages = package_list,
      zip_safe = False,

      install_requires = [
          'astropy>=1.1.2',
          'basemap',
          'boto3>=1.4.4',
          'dispy==4.8.1',
          'dask',
          'graphviz >= 0.7.0',
          'ipython >= 6.0.0',
          'ipywidgets >= 6.0.0',
          'matplotlib >= 2',
          'numpy>=1.10',
          'pandas>=0.17',
          'paramiko >= 2.1.0',
          'psutil >= 5.2',
          'pyproj',
          'pyhht',
          'PyWavelets',
          'scikit-dataaccess>=1.9.0',
          'scikit-learn >= 0.17',
          'scipy',
          'seaborn >= 0.8',
          'setuptools',
          'six >= 1.10.0',
          'statsmodels >= 0.8',
          'tqdm',
      ],

      description = 'A package for Computer-Aided Discovery',
      author = 'MITHAGI',
      author_email='skdaccess@mit.edu',
      classifiers=[
          'Topic :: Scientific/Engineering',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Programming Language :: Python :: 3 :: Only'
          ],
      python_requires='>=3.4',
      url = 'https://github.com/MITHaystack/scikit-discovery',
      package_data={'skdiscovery': ['license/LICENSE','license/LGPL_LICENSE','license/MIT_LICENSE',
                                    'docs/skdiscovery_doxygen.pdf','examples/Amazon_GUI.ipynb',
                                    'examples/Amazon_Offload.ipynb']},

      long_description = readme,
      long_description_content_type='text/markdown'
)
