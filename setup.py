#!/usr/bin/env python
# Time-stamp: <2011-12-29 13:03:08 sunhf>
import sys
import stat
from setuptools import setup, find_packages, Extension


# from distutils.core import setup, Extension,run_setup
from Cython.Distutils import build_ext

# def check_pkg_dependencies():
#     try:
#         import Bio
#     except ImportError, e:
#         sys.stderr.write("CRITICAL: BioPython must be installed\n")
#         sys.exit(1)

def main():
    if not float(sys.version[:3])>=2.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 is recommended!\n")
        sys.exit(1)
    # check_pkg_dependencies()
    ext_modules = [Extension("MPIPE.summary_c",["MPIPE/ext/summary_c.pyx"])]
    setup(name="MPIPE",
          description="A pipeline for motif finding and scoring",
          package_dir={'MPIPE' : 'MPIPE'},
          packages=['MPIPE'],
          install_requires=['biopython', 'rpy2', "cython"],
          scripts=['bin/Mpipe.py','bin/seq_GC_view.py','bin/motif_xml_view.py'],
          author='Hanfei Sun',
          version='0.10',
          author_email='hfsun.tju@gmail.com',
          url='http://samuthing.com',

          classifiers=[
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            ],
          cmdclass = {'build_ext' : build_ext},
          ext_modules = ext_modules
              
          )
 
if __name__ == '__main__':
    main()

