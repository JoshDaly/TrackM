from distutils.core import setup

setup(
    name='TrackM',
    version='0.0.1',
    author='Josh Daly',
    author_email='joshua.daly@uqconnect.edu.au',
    packages=['trackm'],
    scripts=['bin/TrackM'],
    url='http://pypi.python.org/pypi/TrackM/',
    license='GPLv3',
    description='TrackM',
    long_description=open('README.md').read(),
    install_requires=[],
)

