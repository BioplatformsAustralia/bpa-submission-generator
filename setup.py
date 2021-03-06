from setuptools import setup, find_packages

setup(author="CCG, Murdoch University",
      author_email="help@bioplatforms.com",
      description="Generate metadata for international submissions",
      license="GPL3",
      keywords="",
      url="https://github.com/BioplatformsAustralia/bpa-submission-generator",
      name="bpasubmit",
      version="1.2.2",
      packages=find_packages(
          exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
      entry_points={
          'console_scripts': [
              'bpa-submit=bpasubmit.cli:main',
          ],
      })
