from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()
with open("requirements.txt", "r") as f:
    install_requires = f.read().splitlines()

setup(
    name="cartogrampy",
    version="0.0.5",
    description="Python package containing algorithms to calculate various cartogram",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Daniel Lewis, Arturas Eidukas",
    author_email="iwiivi@gmail.com",
    url="https://github.com/danlewis85/cartogrampy",
    install_requires=install_requires,
    packages=find_packages(),
)
