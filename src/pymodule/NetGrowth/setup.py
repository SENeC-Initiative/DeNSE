import os

from setuptools import setup


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as buf:
        return buf.read()


conf = dict(
        name='geometry',
        version='0.1',
        description='geometry module',
        long_description=read('README.md'),
        author='cocconat, tanguy',
        author_email='g33k@paranoici.org',
        license='AGPL',
        packages=['geometry'],
        install_requires=[
            'shapely',
            'descartes',
            'svg.path',
        ],
        zip_safe=False,

        classifiers=[
          "License :: OSI Approved :: GNU Affero General Public License v3",
          "Operating System :: POSIX :: Linux",
          "Programming Language :: Python :: 2",
        ])


if __name__ == '__main__':
    setup(**conf)
