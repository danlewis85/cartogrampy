from setuptools import setup, find_packages
with open('README.md', 'r') as file:
    long_description = file.read()

setup(name='cartogrampy',
      version='0.0.1',
      description='A cartogram plotting package for python.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Daniel Lewis',
      author_email='iwiivi@gmail.com',
      url="https://github.com/danlewis85/cartogrampy",
      packages=find_packages()
    )
