import setuptools

setuptools.setup(
    name="pyamtrack",
    packages=['pyamtrack'],
    package_data={'': ['lib/*']},
    python_requires='>=3.6',
    use_scm_version= {
        "root": "../libamtrack/",
        "relative_to": __file__,
        "local_scheme": "no-local-version"
    },
    setup_requires=['cffi>=1.13.0', 'setuptools_scm'],
    install_requires=['cffi>=1.13.0'],
    author='Leszek Grzanka',
    author_email="leszek.grzanka@ifj.edu.pl",
    description="libamtrack python wrapper",
    url="https://github.com/libamtrack/library",
    long_description = 'libamtrack pythonic wrapper',
    long_description_content_type='text/x-rst',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Physics',

        # OS and env
        'Environment :: Console',
        'Operating System :: POSIX :: Linux',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        # Specify the Python versions you support here.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],
)
