import setuptools

setuptools.setup(
    name="pyamtrack",
    version="0.1.0",
    packages=['pyamtrack'],
    package_data={'': ['lib/*']},
    python_requires='>=3.6',
    author_email="leszek.grzanka@ifj.edu.pl",
    description="libamtrack python wrapper",
    url="https://github.com/libamtrack/library")
