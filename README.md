![cartogrampy logo](https://i.imgur.com/W1XfhLB.png "cartogrampy logo")
# Cartogrampy 
Cartogram and related algorithms for GeoPandas

## Installation

### Using Pip

Prerequisited 
- python3
- git
- make

To install this package using pip run these commands:
``` sh
git clone git@github.com:danlewis85/cartogrampy.git
cd cartogrampy
make install
pip install . 
```

If you know what you are doing and/or want to use `pipenv`,
then run this sequence of commands:

``` sh
git clone git@github.com:danlewis85/cartogrampy.git
cd cartogrampy
make pipenv
pipenv install . 
```

If you are unfortunate enough to exist in Windows,
run this sequence of commands in your command prompt:

``` sh
git clone git@github.com:danlewis85/cartogrampy.git
cd cartogrampy
pip install --upgrade pipenv 
pipenv install -r requirements.txt
pipenv install .
```
