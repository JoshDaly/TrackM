from distutils.core import setup

setup(
    name='TrackM',
    version='0.2.0',
    author='Josh Daly, Michael Imelfort',
    author_email='joshua.daly@uqconnect.edu.au',
    packages=['trackm'],
    scripts=['bin/trackm'],
    url='http://pypi.python.org/pypi/TrackM/',
    license='GPLv3',
    description='TrackM - track HGT in microbial genomes',
    long_description=open('README.md').read(),
    install_requires=[
        "jsonpickle >= 0.7.1",
        "DancingPeasant >= 0.1.3"
        "pyzmq >= 14.3.1"
    ],
)

