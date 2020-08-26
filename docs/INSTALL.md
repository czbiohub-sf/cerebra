Installation Instructions
============

Dependencies
------------
`cerebra` depends on some (fairly standard) packages and libraries. 
Before installing it might be a good idea to make sure all of the requisite packages are installed on your system (_note:_ if installing with Docker you can skip this step). 

__MacOS Dependencies:__
```
sudo pip install setuptools
brew update
brew install openssl
brew install zlib
```

__Linux Dependencies:__
```
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

As of present `cerebra` is not installable on Windows. 
`cerebra` depends on the [`pysam`](https://pysam.readthedocs.io/en/latest/index.html) library (or rather, `pysam` is a dependency-of-a-dependency) and currently this library is only available on Unix-like systems. 
Windows solutions like [WSL](https://docs.microsoft.com/en-us/windows/wsl/) exist for overcoming precisely this challange, however, `cerebra` has not been tested on WSL or any other Unix-like subsystem for Windows.    


Installation (for users)
------------
There are four different methods available to install `cerebra`.
Choose one of the following:

__With [Docker](https://hub.docker.com/r/lincolnharris/cerebra) (recommended)__          
```
docker pull lincolnharris/cerebra
docker run -it lincolnharris/cerebra
```      
_warning: this image will take up ~1Gb of space._               

__With traditional git clone and the python standard library [`venv`](https://docs.python.org/3/library/venv.html) module__
```
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
python3 -m venv cerebra-dev
source cerebra-dev/bin/activate
pip3 install [--user] . 
```

__With traditional git clone and [conda](https://docs.conda.io/en/latest/)__
``` 
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install [--user] . 
```

__From [PyPi](https://pypi.org/project/cerebra/)(system-wide installation, NOT RECOMMENDED)__           
For novice users, its generally a better idea to install packages within virtual environments. 
However, `cerebra` can be installed system-wide, if you so choose. 
```
pip install cerebra

# OR, if you dont have root privileges
pip install --user cerebra
```


Installation (for developers)
------------
Here's how to set up cerebra for local development. 
After installing the requisite [dependencies](https://github.com/czbiohub/cerebra/blob/master/docs/INSTALL.md#dependencies):

1.  Fork the `cerebra` repo on GitHub: https://github.com/czbiohub/cerebra
2.  Clone your fork locally:

        $ git clone https://github.com/your-name/cerebra.git

3.  Install your local copy into a virtualenv. Using the standard library [`venv`](https://docs.python.org/3/library/venv.html) module: 

        $ cd cerebra
        $ python3 -m venv cerebra-dev
        $ source cerebra-dev/bin/activate
        $ pip install -r requirements.txt -r test_requirements.txt -e . 

4.  Create a branch for local development:

        $ git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

5.  When you're done making changes, check that your changes pass flake8 and the tests:

        $ make test
        $ make coverage
        $ make lint

6.  Commit your changes and push your branch to GitHub:

        $ git add .
        $ git commit -m "Your detailed description of your changes."
        $ git push origin name-of-your-bugfix-or-feature

7.  Submit a pull request through the GitHub website.
See [CONTRIBUTING.md](https://github.com/czbiohub/cerebra/blob/master/docs/CONTRIBUTING.md) for more. 
