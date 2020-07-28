Installation
------------
There are four different methods available to install `cerebra`.

__With [Docker](https://hub.docker.com/r/lincolnharris/cerebra)__                
```
docker pull lincolnharris/cerebra
```                 

__With traditional git clone and the python standard library [venv](https://docs.python.org/3/library/venv.html) module__
```
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
python3 -m venv cerebra-dev
source cerebra-dev/bin/activate
pip3 install -e . 
```

__With traditional git clone and [conda](https://docs.conda.io/en/latest/)__
``` 
git clone https://github.com/czbiohub/cerebra.git
cd cerebra
conda create -n cerebra python=3.7
conda activate cerebra
pip3 install -e . 
```

__From [PyPi](https://pypi.org/project/cerebra/) (system-wide installation)__              
```
pip install cerebra

# OR, if you dont have root privileges
pip install --user cerebra
```

`cerebra` depends on some (fairly standard) packages and libraries. If youre having trouble installing, it might be a good idea to make sure you have all of the requisite dependendices installed first (_note:_ if installing with Docker you can skip this step). 

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

As of present `cerebra` is not installable on Windows. `cerebra` depends on the [`pysam`](https://pysam.readthedocs.io/en/latest/index.html) library -- or rather, `pysam` is a dependency-of-a-dependency -- and currently this library is only available on Unix-like systems. 
