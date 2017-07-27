from setuptools import setup, find_packages

package_name = 'skdiscovery'

package_list = find_packages()
package_list = [ package_name + '.' + package for package in package_list]
package_list.insert(0, package_name)


setup(name     = package_name,
      version  = '1.0',
      packages = package_list,
      package_dir = {package_name: ''},
      zip_safe = False,
      install_requires = ['setuptools',
                          'numpy',
                          'scipy',
                          'pandas',
                          'basemap',
                          'sklearn',
                          'skdaccess']
       )
 
