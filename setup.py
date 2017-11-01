from setuptools import setup
from setuptools import find_packages

package_name = 'scikit-discovery'

package_list = find_packages()

setup(name     = package_name,
      version  = '0.9.8',
      packages = package_list,
      zip_safe = False,
      
      install_requires = ['basemap',
                          'tqdm',
                          'numpy>=1.10',
                          'pandas>=0.17',
                          'scipy',
                          'setuptools',
                          'astropy>=1.1.2',
                          'scikit-dataaccess>=1.9.0',
                          'psutil>=5',
                          'boto3>=1.4.4',
                          'statsmodels >= 0.8',
                          'graphviz >= 0.7.0',
                          'paramiko >= 2.1.0',
                          'matplotlib >= 2',
                          'dispy==4.8.1',
                          'ipython >= 6.0.0',
                          'ipywidgets >= 6.0.0',
                          'psutil >= 5.2',
                          'seaborn >= 0.8',
                          'six >= 1.10.0',
                          'scikit_learn >= 0.17',
                          'PyWavelets',
                          'pyhht'],
      
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
      package_data={'skdiscovery': ['license/LICENSE','license/LGPL_LICENSE','license/MIT_LICENSE']},

)
