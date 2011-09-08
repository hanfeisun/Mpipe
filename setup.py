import sys
import stat
from distutils.core import setup, Extension
def check_pkg_dependencies():
    print "lala"
    try:
        import Bio
    except ImportError, e:
        sys.stderr.write("CRITICAL: BioPython must be installed\n")
        sys.exit(1)

def main():
    print "ok"
    if not float(sys.version[:3])>=2.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 is recommended!\n")
        sys.exit(1)
    # check_pkg_dependencies()
    setup(name="Forscholoarshipexam",
          description="For scholarship exam",
          package_dir={'fse' : 'lib'},
          packages=['fse'],
          scripts=['bin/pipe.py','bin/seq_GC_view.py','bin/motif_xml_view.py'],
          author='Hanfei Sun',
          author_email='hfsun.tju@gmail.com',
          url='http://samuthing.com',

          classifiers=[

            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            ]
          )
 
if __name__ == '__main__':
    main()

